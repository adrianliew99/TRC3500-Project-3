"""
plot_fusion.py - Event-level breath-cycle fusion of S1 (strain) and S2 (thermistor).

On launch the script enters a CALIBRATION PHASE that prints S1's raw
ADC voltage live to the console. The user adjusts the strain-band
tension while watching the bridge output settle, then presses any
key (Windows) to start the live plot, which begins with a 10 s
filter warmup (controlled by WARMUP_S) during which no breaths are
detected.

Reads BOTH ADC channels from each 'DATA|S1:<int>|S2:<int>' UART line
and runs them through independent processing chains:

  1. Long rolling-average drift removal.
  2. Short rolling-average noise smoothing.
  3. Hilbert-phase peak detector. A sliding HILBERT_BUF_S-second
     window of the bandpass signal is run through an FFT-based
     analytic-signal transform every sample. A breath event is
     emitted when the instantaneous phase wraps from + to - through
     zero (i.e. the bandpass is at its peak) at the lookback index,
     gated by an envelope floor (the adaptive threshold still shown
     on the panels) and a MIN_BREATH_S refractory.

The Hilbert-based detector replaces the older above-/below-threshold
state machine. Phase wraps exactly once per breath even when the
waveform is asymmetric or carries sub-cycle motion ripple, so
detection no longer doubles up at slow rates the way a naive peak
threshold did. See `analyze_dsp.py` for the per-method comparison
that motivated the switch.

The two sensors are out of phase and have COMPLEMENTARY failure
modes, so adding their waveforms makes no sense. Instead we pick
ONE sensor at a time as the authoritative breath source, with a
THREE-step mode ladder driven by the median of recent S1-to-S1
intervals (S1 fires per-breath in every mode, so its interval
history is the cleanest rate gauge):

SLOW mode (median S1 interval >= MODERATE_MODE_INTERVAL_S):
  Trust S2 (thermistor) 100 %. Each S2 trough is taken directly as
  one verified breath. The TMP61 bead has enough time to follow
  each breath fully. Both thresholds run at their nominal floors.

MODERATE mode (FAST_MODE_INTERVAL_S <= median < MODERATE_MODE_INTERVAL_S):
  Still trust S2, but breaths are smaller and the thermistor
  swing per breath shrinks. S2's adaptive threshold is scaled down
  by MODERATE_MODE_S2_THR_SCALE so the shallower troughs still
  cross. S1 keeps its nominal threshold (it isn't sourcing breaths
  in this mode, just driving the mode decision).

FAST mode (median S1 interval < FAST_MODE_INTERVAL_S):
  Trust S1 (strain) 100 %. Each S1 peak is one verified breath. S2
  can't keep up at this rate (the bead's thermal mass averages
  multiple breaths). S1's adaptive threshold is scaled down by
  FAST_MODE_S1_THR_SCALE for the same reason - faster breathing
  tends to be shallower.

Mode switching is automatic and per-sample.

A breathing CYCLE is the span between two consecutive verified
breaths, drawn as a shaded region in panel 3.

Figure layout (live, scrolling WINDOW_S seconds):
  Panel 1: filtered S1 (strain) bandpass + threshold + per-sensor breath markers
  Panel 2: filtered S2 (thermistor) bandpass + threshold + per-sensor breath markers
  Panel 3: verified-breath timeline. Faint normalised S1/S2 traces
           sit underneath as context. Each green vertical line is a
           SLOW-mode breath sourced from S2; each orange vertical
           line is a FAST-mode breath sourced from S1. Translucent
           green regions bracket each cycle - the span from one
           verified breath to the next.

Firmware contract: 'DATA|S1:<int>|S2:<int>\\n' lines at ~100 Hz, 460800 baud.
"""

import sys
import time
from collections import deque

import numpy as np
import serial
import serial.tools.list_ports
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection


def _hilbert_analytic(x_arr):
    """Return the analytic signal (real + j*Hilbert) of a real 1-D array via FFT.

    The instantaneous phase angle(z) wraps exactly ONCE per cycle of the
    dominant frequency in x_arr, even when the underlying waveform has
    harmonic distortion (asymmetric inhale/exhale, motion ripple on the
    strain band, etc). That is the property the per-sensor detector below
    exploits: it ignores sub-cycle wobbles that a peak/threshold detector
    would otherwise count as extra breaths.
    """
    n = len(x_arr)
    X = np.fft.fft(x_arr)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1:(n + 1) // 2] = 2.0
    return np.fft.ifft(X * h)

try:
    import msvcrt           # Windows non-blocking keyboard, used by the
                            # pre-plot calibration phase. The rest of the
                            # script is platform-agnostic.
except ImportError:
    msvcrt = None

try:
    import winsound         # Windows beep used to mark the moment the
                            # live plot starts streaming, so the user
                            # has an audible cue that recording is on.
except ImportError:
    winsound = None


# --- Link to firmware -----------------------------------------------------
BAUD_RATE        = 460800
EMIT_HZ_NOMINAL  = 100

# --- Plot ----------------------------------------------------------------
WINDOW_S         = 100
UPDATE_EVERY     = 10

# --- Per-sensor filter + detector knobs (taken from the single-channel
#     scripts so behaviour matches what those plots already show) --------
# S1 = strain band on the belly. Polarity not inverted.
S1_FILTER_WIN_S  = 5
S1_LP_WIN_S      = 0.3
S1_FLOOR         = 150         # counts (~7.5 mV)
S1_INVERT        = False
# S2 = thermistor near the face. Cooling -> ADC drop, so invert so
# inhales appear as positive peaks (same as plot_adc_s2.py).
S2_FILTER_WIN_S  =5
S2_LP_WIN_S      = 0.5
S2_FLOOR         = 100         # counts (~5 mV)
S2_INVERT        = True

# Shared peak-detector knobs.
# Envelope/threshold parameters still drive the +/-thr lines drawn on each
# per-sensor panel AND act as the amplitude floor that gates Hilbert-phase
# events (so a quiet baseline cannot fire spurious breaths). They no longer
# drive the event-detection state machine itself - that is now Hilbert.
ENVELOPE_TAU_S   = 1.5
ADAPTIVE_K       = 0.12
MIN_BREATH_S     = 0.8
# Longer refractory used by the parallel simple-threshold detector so
# it doesn't fire on slow-rate sub-cycle motion ripple. 1.0 s admits
# any breath up to 60 BPM and lets the simple detector catch the FIRST
# post-warmup peak even if it lands right at the warmup boundary.
SIMPLE_REFRACTORY_S = 1.0
# With calibration-baseline prefill the HP/LP rolling means are at
# steady state in one sample, so the only thing still warming up is
# the Hilbert sliding buffer (HILBERT_BUF_S). Set warmup just longer
# than that so the buffer is full of real samples by the time
# detection is allowed.
WARMUP_S         = 20.0

# --- Hilbert-phase detector ---------------------------------------------
# Each per-sensor detector keeps a sliding buffer of bandpass samples and,
# every sample after the buffer fills, computes the analytic signal via
# FFT. A breath event is emitted when the instantaneous phase wraps from
# positive to negative (= bandpass peak) AND the analytic envelope at
# that moment exceeds the adaptive threshold AND the per-sensor refractory
# (MIN_BREATH_S) has elapsed.
#
# HILBERT_BUF_S      length of the sliding analysis window. Must cover at
#                    least one full breath at the slowest expected rate
#                    so the FFT can resolve the breathing frequency.
#                    Tuned to 5 s: that is exactly one cycle at the
#                    slowest target speed (12 BPM metronome => 5 s breath
#                    period) AND about 3.7 cycles at the fastest target
#                    speed (80 BPM metronome => 1.36 s breath period).
#                    Originally 8 s, but at that length a single startup
#                    transient on S1 (the huge spike when the rolling
#                    mean first sees a real breath) dominated the FFT for
#                    eight seconds afterwards, pinning the apparent
#                    dominant frequency low and stalling phase rotation
#                    so subsequent fast breaths could not fire. Shorter
#                    buffer => recent breaths dominate the spectrum.
#                    Trade-off: cycles longer than 5 s (e.g. super_slow
#                    at ~7.5 s) are no longer resolvable; that mode is
#                    not in the test target list per project memory.
# HILBERT_LOOKBACK_S phase is read this far back from the buffer's leading
#                    edge so we avoid the FFT's edge artefacts. Sets the
#                    event-reporting latency.
HILBERT_BUF_S      = 3.0    # was 5.0; shrunk so a fresh breath dominates
                            # the FFT spectrum within ~one cycle instead of
                            # being averaged out by older samples
HILBERT_LOOKBACK_S = 0.15   # half what it was; less event-emit latency

# --- Fusion --------------------------------------------------------------
# Three-mode ladder driven by the median of recent S1-to-S1 intervals.
# SLOW     (interval >= MODERATE_MODE_INTERVAL_S):
#   S2 sources every breath, both thresholds at nominal.
# MODERATE (FAST_MODE_INTERVAL_S <= interval < MODERATE_MODE_INTERVAL_S):
#   S2 still sources, but its threshold is scaled down because
#   the per-breath swing has shrunk with the faster rate.
# FAST     (interval < FAST_MODE_INTERVAL_S):
#   S1 sources every breath; the thermistor can't follow individual
#   breaths this fast. S1's threshold is scaled down too.
# Mode switching is automatic and per-sample. A breathing CYCLE is
# the span between two consecutive verified breaths (shaded region
# in panel 3).
# Mode boundaries on the median of recent S1-to-S1 intervals:
#     median  >=  MODERATE_MODE_INTERVAL_S  ->  SLOW
#     FAST_MODE_INTERVAL_S  <=  median  <  MODERATE_MODE_INTERVAL_S -> MODERATE
#     median  <  FAST_MODE_INTERVAL_S       ->  FAST
MODERATE_MODE_INTERVAL_S = 3.0
FAST_MODE_INTERVAL_S     = 1.1
FAST_MODE_LOOKBACK       = 3    # number of recent S1 intervals to median
# Threshold scaling per mode. When breathing is faster, the per-breath
# excursion in the trusted sensor is smaller, so we shrink that
# sensor's adaptive threshold to keep catching the crossings.
# Lower = more sensitive (more breaths caught, more noise risk).
#   SLOW     : S1 scale 1.0, S2 scale 1.0
#   MODERATE : S1 scale 1.0, S2 scale MODERATE_MODE_S2_THR_SCALE  (S2 trusted)
#   FAST     : S1 scale FAST_MODE_S1_THR_SCALE, S2 scale 1.0       (S1 trusted)
MODERATE_MODE_S2_THR_SCALE = 0.5
FAST_MODE_S1_THR_SCALE     = 0.5

# --- ADC conversion ------------------------------------------------------
ADC_FULLSCALE    = 65535
VREF_V           = 3.3


def calibration_phase(ser):
    """Stream S1 raw voltage to the console until the user presses a key.

    Lets the user adjust the strain-band tension while watching how
    the Wheatstone-bridge output settles. A reasonable target is a
    stable reading somewhere near mid-rail (~1.65 V on a 3.3 V
    reference) so each breath has headroom to swing in both
    directions before clipping the bridge.

    Returns when any key is pressed (Windows: msvcrt). On other
    platforms keypress detection isn't wired up, so the loop runs
    until Ctrl+C - the rest of the script still works.

    Returns (baseline_s1, baseline_s2): mean of the LAST few seconds of
    samples observed during calibration. These get passed into the
    SensorChain prefill so the HP rolling-mean starts at the real
    breathing baseline rather than the first arbitrary sample - which
    eliminates the slow DC drift that otherwise lingers for ~5 s
    after warmup ends and causes the first breath after warmup to be
    missed.
    """
    print()
    print("--- CALIBRATION PHASE ---")
    print("Adjust the strain-band tension.")
    print()
    if msvcrt is None:
        print("[!] Non-Windows: keypress detect unavailable. Ctrl+C to abort.")
    else:
        print("Press any key to start the live plot (10 s warmup applies).")
    print()

    # Keep a fixed-length buffer of the most recent (y1, y2) samples
    # observed during calibration. On exit we average it to get a
    # baseline estimate that's a much better prefill value for the HP
    # rolling mean than just "the first sample we ever saw".
    baseline_n = int(5.0 * EMIT_HZ_NOMINAL)   # last ~5 s of samples
    y1_hist = deque(maxlen=baseline_n)
    y2_hist = deque(maxlen=baseline_n)

    last_print_t = 0.0
    while True:
        if msvcrt is not None and msvcrt.kbhit():
            msvcrt.getch()        # consume the key
            print()               # leave the live-line behind
            print("[+] Calibration confirmed. Starting live plot...\n")
            bs1 = (sum(y1_hist) / len(y1_hist)) if y1_hist else None
            bs2 = (sum(y2_hist) / len(y2_hist)) if y2_hist else None
            return bs1, bs2

        line = ser.readline().decode("utf-8", errors="ignore").strip()
        if not line.startswith("DATA"):
            continue
        fields = {kv.split(":")[0]: kv.split(":")[1]
                  for kv in line.split("|")[1:] if ":" in kv}
        if "S1" not in fields:
            continue
        try:
            y1 = int(fields["S1"])
            y2 = int(fields["S2"]) if "S2" in fields else None
        except ValueError:
            continue

        y1_hist.append(y1)
        if y2 is not None:
            y2_hist.append(y2)

        v = y1 * VREF_V / ADC_FULLSCALE

        # Throttle the carriage-return live update to ~10 Hz so the
        # terminal stays readable and the kbhit() check still runs
        # often. ser.readline() at 100 Hz drives the outer rate.
        now = time.time()
        if now - last_print_t > 0.1:
            last_print_t = now
            print(f"  S1 raw: {y1:6d} counts  ({v:5.3f} V)  "
                  f"[press any key to continue]   ",
                  end="\r", flush=True)


def auto_detect_com_port():
    # Same idiom used by the single-channel scripts.
    print("--- Auto-Detecting STM32 ---")
    for device in serial.tools.list_ports.comports():
        desc = str(device)
        if "STMicroelectronics" in desc or "STLink" in desc:
            port = desc.split()[0]
            print(f"[+] STM detected on port: {port}")
            return port
    print("[-] Could not detect STM32 - check USB cable and ST-Link drivers.")
    sys.exit(1)


class SensorChain:
    """Per-channel filter + Hilbert-phase cycle detector.

    Streams ADC samples through:
      - rolling-average HP (drift removal)
      - rolling-average LP (HF noise smoothing)
      - sliding-window FFT Hilbert -> instantaneous phase + envelope

    Emits a breath event (peak_t, end_t) every time the analytic
    phase wraps from + to - through zero, gated by the adaptive
    envelope floor and MIN_BREATH_S refractory. peak_t == end_t for
    this detector; end_t is retained only for interface symmetry
    with the fuser, which has historically taken both. See the
    module docstring for the motivation behind the Hilbert switch.
    """

    def __init__(self, name, fs_hz, filter_win_s, lp_win_s, floor, invert,
                 baseline_y=None):
        self.name    = name
        self.invert  = invert
        self.floor   = floor
        self.baseline_y = baseline_y    # used for the HP-mean prefill

        # Long rolling-average ("highpass") to subtract slow drift.
        self.hp_n    = int(filter_win_s * fs_hz)
        self.hp_win  = deque(maxlen=self.hp_n)
        self.hp_sum  = 0.0

        # Short rolling-average ("lowpass") to round off fast noise.
        self.lp_n    = max(2, int(lp_win_s * fs_hz))
        self.lp_win  = deque(maxlen=self.lp_n)
        self.lp_sum  = 0.0

        # Envelope follower the adaptive threshold rides on top of.
        # Also used to envelope-normalise the channel before drawing
        # it in the fused panel.
        self.env_decay = 1.0 - 1.0 / (ENVELOPE_TAU_S * fs_hz)
        self.envelope  = 0.0
        self.threshold = floor
        # Multiplier applied to the final threshold each sample. The
        # main loop flips this between 1.0 (normal) and
        # FAST_MODE_S1_THR_SCALE so the fuser's mode decision can
        # boost detector sensitivity when fast mode kicks in.
        self.threshold_scale = 1.0

        # Per-sensor Hilbert-phase detector state.
        # Sliding buffers of the most recent bandpass samples and their
        # timestamps. When full, every new sample triggers an FFT-based
        # analytic-signal computation; a phase wrap through zero at the
        # lookback index marks one breath peak. Deliberately NO window
        # is applied before the FFT: a Hamming or Hann window attenuates
        # the buffer edges by 90+ %, so the phase reading at the lookback
        # index (which lives near the edge) would be dominated by noise.
        # Spectral-leakage artefacts from the rectangular window are
        # tolerated; the dominant breathing-frequency phase is preserved.
        self.hilb_buf_n = max(8, int(HILBERT_BUF_S * fs_hz))
        self.hilb_lk_n  = max(1, int(HILBERT_LOOKBACK_S * fs_hz))
        self.bp_history = deque(maxlen=self.hilb_buf_n)
        self.t_history  = deque(maxlen=self.hilb_buf_n)
        self.prev_phase = None
        # Parallel simple-threshold detector state. Hilbert needs a clean
        # periodic signal to lock onto; at super_fast rates the strain
        # band's per-breath swing is small enough that the spectrum lacks
        # a clear dominant peak, and Hilbert misses cycles. The simple
        # detector below complements Hilbert by firing on the
        # bp-upward-crosses-floor moment, which is what the original
        # plot_fusion threshold detector did. Both detectors share the
        # MIN_BREATH_S refractory so they cannot double-count a breath.
        self.simple_above    = False
        self._simple_peak_bp = 0.0
        # Refractory + per-sensor peak history (used for the red dashed
        # markers in the per-sensor panel AND for the fuser's interval gauge).
        self.peak_time   = 0.0
        self.last_peak_t = -1e9
        self.breath_times = []

    def step(self, y, now_t):
        """Process one ADC sample.
        Returns (bp_in_counts, threshold_in_counts, cycle_event_or_None).
        cycle_event = (start_t, end_t) when an inhale just completed."""
        # Pre-fill the rolling-mean buffers on the first call. Without
        # this the HP rolling mean grows from 1 to hp_n samples, so for
        # the first several seconds the mean lags the true baseline and
        # drags bp off zero - a slow DC drift that swallows the first
        # post-warmup breath. We prefill with `baseline_y` (the average
        # observed during calibration, when the user is breathing
        # normally) if available, falling back to the very first live
        # sample otherwise. The calibration mean is the better choice:
        # the FIRST sample might land mid-inhale or mid-exhale and bias
        # the mean by a full breath amplitude, whereas the calibration
        # mean already averages over many breath cycles.
        if not self.hp_win:
            seed = float(self.baseline_y) if self.baseline_y is not None else float(y)
            self.hp_win.extend([seed] * self.hp_n)
            self.hp_sum = seed * self.hp_n
            self.lp_win.extend([0.0] * self.lp_n)
            self.lp_sum = 0.0

        self.hp_sum -= self.hp_win[0]
        self.hp_win.append(y); self.hp_sum += y
        hp_val = y - self.hp_sum / len(self.hp_win)

        self.lp_sum -= self.lp_win[0]
        self.lp_win.append(hp_val); self.lp_sum += hp_val
        bp = self.lp_sum / len(self.lp_win)
        if self.invert:
            bp = -bp

        self.envelope  = max(abs(bp), self.envelope * self.env_decay)
        self.threshold = self.threshold_scale * max(self.floor,
                                                    ADAPTIVE_K * self.envelope)

        # Push the new bandpass sample + timestamp into the sliding
        # Hilbert window. The detector only fires once the window is
        # full so the FFT has enough samples to resolve the breathing
        # frequency; before then we are still in warmup.
        self.bp_history.append(bp)
        self.t_history.append(now_t)

        event = None
        if now_t >= WARMUP_S and len(self.bp_history) == self.hilb_buf_n:
            arr = np.asarray(self.bp_history, dtype=float)
            z = _hilbert_analytic(arr)
            # Read phase/envelope hilb_lk_n samples back from the leading
            # edge - the FFT distorts the very edge of the window, but
            # has settled by the lookback point. Cost: each detected event
            # is reported HILBERT_LOOKBACK_S after it actually happened,
            # which is fine because it shifts all events uniformly so the
            # fuser's interval gauge is unaffected.
            idx  = self.hilb_buf_n - 1 - self.hilb_lk_n
            phi  = float(np.angle(z[idx]))
            env  = float(np.abs(z[idx]))
            t_at = self.t_history[idx]

            # Peak = phase crosses zero UPWARD. Convention reminder: for a
            # real cosine x(t)=cos(omega*t), the analytic signal is
            # exp(j*omega*t), whose phase INCREASES linearly with t. So
            # phase=0 corresponds to the bandpass peak, and the cycle
            # advances from negative phase through zero to positive phase.
            # The small-step guard rules out the spurious -pi <-> +pi wrap
            # (a trough). MIN_BREATH_S refractory caps the event rate at
            # 60/0.8 = 75 BPM. The amplitude floor uses the analytic
            # envelope vs the bare-counts floor (not the bp-peak-driven
            # adaptive threshold, which a single big breath would inflate
            # and then take seconds to decay, locking out smaller
            # subsequent breaths). The +/- threshold lines on the panels
            # still come from the bp-peak envelope so the user sees the
            # same visual cue they always did.
            if (self.prev_phase is not None
                    and self.prev_phase < 0.0 and phi >= 0.0
                    and (phi - self.prev_phase) < np.pi
                    and env > (self.floor * self.threshold_scale)
                    and (t_at - self.last_peak_t) >= MIN_BREATH_S):
                self.last_peak_t = t_at
                self.peak_time   = t_at
                self.breath_times.append(t_at)
                # peak_t == end_t for the Hilbert detector: every event is
                # anchored at the phase=0 crossing. The fuser only uses
                # peak_t; end_t is retained for interface compatibility.
                event = (t_at, t_at)
            self.prev_phase = phi

        # Parallel simple-threshold detector. Runs unconditionally; it
        # cannot be gated on fuser mode because reaching fast mode
        # requires S1 events, but at super_fast the Hilbert detector
        # silently fails to produce them - chicken and egg. Uses its own
        # longer refractory (SIMPLE_REFRACTORY_S = 1.2 s) so it does NOT
        # fire on slow-rate sub-cycle motion ripple (which sits ~0.8-1 s
        # apart and slips past the standard MIN_BREATH_S=0.8 s gate).
        # 1.2 s still admits any breath up to 50 BPM - the band of rates
        # where Hilbert is the one that struggles. Shared last_peak_t
        # with Hilbert prevents the two paths from double-counting.
        if (now_t >= WARMUP_S and event is None):
            simple_thr = self.floor * self.threshold_scale
            if not self.simple_above:
                if bp > simple_thr:
                    self.simple_above = True
                    self.peak_time   = now_t
                    self._simple_peak_bp = bp
            else:
                if bp > self._simple_peak_bp:
                    self._simple_peak_bp = bp
                    self.peak_time       = now_t
                if bp < simple_thr:
                    self.simple_above = False
                    if (self.peak_time - self.last_peak_t) >= SIMPLE_REFRACTORY_S:
                        self.last_peak_t = self.peak_time
                        self.breath_times.append(self.peak_time)
                        event = (self.peak_time, now_t)

        return bp, self.threshold, event


class FusedEventDetector:
    """Mode-switching breath-event detector.

    The two sensors have complementary strengths:
      - S2 (thermistor) is sluggish but produces a clean trough on
        each breath when it has enough time to respond. It is the
        better source when breathing is slow.
      - S1 (strain) is instant but noisier (motion artefacts pick up
        as mistriggers). It is the better source when breathing is
        fast enough that S2 cannot keep up.

    So we pick ONE sensor at a time as the authoritative source for
    breath events, switching between them based on the median of the
    most recent S1-to-S1 intervals (S1 always fires per-breath in
    either mode, so its interval history is the cleanest gauge of
    the current breathing rate):

      SLOW mode (median S1 interval >= FAST_MODE_INTERVAL_S):
        Each S2 fire is emitted directly as a verified breath. S1
        fires are ignored for emission but still update the interval
        history (so we can detect when breathing speeds up).

      FAST mode (median S1 interval <  FAST_MODE_INTERVAL_S):
        Each S1 fire is emitted directly as a verified breath. S2
        fires are ignored.

    Mode transitions are automatic and silent. Startup defaults to
    SLOW mode until enough S1 fires have accumulated to compute a
    median.

    A breathing CYCLE is the span between two consecutive verified
    breaths (drawn as a shaded region in panel 3).

    Returns an event tuple
        (idx, peak_t, end_t, interval_or_None, bpm_or_None, contrib)
    from ingest(), or None. tick() always returns None (kept for
    interface symmetry; this detector has no time-based state).
    """

    def __init__(self):
        # S1-to-S1 interval history (drives mode selection).
        self.last_s1_t        = None
        self.recent_intervals = deque(maxlen=FAST_MODE_LOOKBACK)
        # Last VERIFIED breath time, for BPM computation.
        self.last_verified_t  = None
        self.cycle_idx        = 0
        # Verified-breath history kept for panel-3 markers + spans.
        self.cycles_start     = []   # verified peak times
        self.cycles_end       = []   # source-event end_t (for span end)
        self.cycles_paired    = []   # True if S2-sourced, False if S1-sourced

    def current_mode(self):
        """Return 'slow', 'moderate', or 'fast' based on median S1 interval.

        - 'slow'     : interval >= MODERATE_MODE_INTERVAL_S  (S2 sourced, full thresholds)
        - 'moderate' : FAST_MODE_INTERVAL_S <= interval < MODERATE_MODE_INTERVAL_S
                       (still S2 sourced, but S2 threshold scaled down)
        - 'fast'     : interval < FAST_MODE_INTERVAL_S       (S1 sourced, S1 threshold scaled down)
        Defaults to 'slow' until enough S1 intervals exist for a median.
        """
        if len(self.recent_intervals) < 2:
            return "slow"
        s = sorted(self.recent_intervals)
        median = s[len(s) // 2]
        if median < FAST_MODE_INTERVAL_S:
            return "fast"
        if median < MODERATE_MODE_INTERVAL_S:
            return "moderate"
        return "slow"

    def _emit(self, peak_t, end_t, sourced_from_s2):
        interval = None
        bpm      = None
        if self.last_verified_t is not None:
            interval = peak_t - self.last_verified_t
            if interval > 0:
                bpm = 60.0 / interval
        self.last_verified_t = peak_t
        self.cycle_idx += 1
        self.cycles_start.append(peak_t)
        self.cycles_end.append(end_t)
        self.cycles_paired.append(sourced_from_s2)
        contrib = ["S2"] if sourced_from_s2 else ["S1"]
        return (self.cycle_idx, peak_t, end_t, interval, bpm, contrib)

    def ingest(self, sensor_id, peak_t, end_t):
        if sensor_id == "S1":
            # Always update S1 interval history, even in slow mode,
            # so we can notice when the user starts breathing faster
            # and switch into fast mode on the spot.
            if self.last_s1_t is not None:
                self.recent_intervals.append(peak_t - self.last_s1_t)
            self.last_s1_t = peak_t
            if self.current_mode() == "fast":
                return self._emit(peak_t, end_t, sourced_from_s2=False)
            return None     # slow / moderate: trust S2, ignore S1 for emission

        # sensor_id == "S2"
        if self.current_mode() == "fast":
            return None     # fast mode: trust S1, ignore S2 for emission
        return self._emit(peak_t, end_t, sourced_from_s2=True)



def main():
    port = auto_detect_com_port()
    try:
        ser = serial.Serial(port, BAUD_RATE, timeout=1)
    except Exception as e:
        print(f"[-] Failed to open {port}: {e}")
        sys.exit(1)
    print(f"[+] Opened {port} @ {BAUD_RATE} baud")

    # Same boot-up handshake as the single-channel scripts: flush any
    # stale bytes and wait briefly for the firmware to start streaming.
    ser.reset_input_buffer()
    deadline = time.time() + 3.0
    saw_data = False
    while time.time() < deadline:
        line = ser.readline().decode("utf-8", errors="ignore").strip()
        if not line:
            continue
        saw_data = True
        if line.startswith("DATA") or line.startswith("STM32_READY"):
            break
    if not saw_data:
        print("[!] No data. Press the BLACK RESET button on the Nucleo.")

    # Pre-plot calibration: lets the user dial in strain-band tension
    # by watching S1 raw voltage in the terminal. Any keypress
    # advances to the live plot.
    baseline_s1, baseline_s2 = calibration_phase(ser)
    # Flush any data accumulated during calibration so the live plot's
    # first samples and the 10 s warmup begin from a fresh buffer.
    ser.reset_input_buffer()

    s1 = SensorChain("S1", EMIT_HZ_NOMINAL, S1_FILTER_WIN_S, S1_LP_WIN_S, S1_FLOOR, S1_INVERT,
                     baseline_y=baseline_s1)
    s2 = SensorChain("S2", EMIT_HZ_NOMINAL, S2_FILTER_WIN_S, S2_LP_WIN_S, S2_FLOOR, S2_INVERT,
                     baseline_y=baseline_s2)
    fuser = FusedEventDetector()

    # --- Build the live figure ------------------------------------------
    plt.ion()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 9), sharex=True)

    line_s1, = ax1.plot([], [], color="tab:green", linewidth=0.8)
    thr1_p,  = ax1.plot([], [], color="tab:orange", linestyle=":", linewidth=0.9, alpha=0.7, label="+/-thr")
    thr1_n,  = ax1.plot([], [], color="tab:orange", linestyle=":", linewidth=0.9, alpha=0.7)
    ax1.axhline(0, color="black", linewidth=0.5, alpha=0.4)
    ax1.set_ylabel("S1 bandpass (mV)")
    ax1.set_title(f"Event-level fusion - live (last {WINDOW_S} s)"
                  f"  |  green = S1 (peak = end inhale)"
                  f"  |  purple = S2 raw (trough = end inhale)"
                  f"  |  red dashed = end of inhale")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right", fontsize=8)

    line_s2, = ax2.plot([], [], color="tab:purple", linewidth=0.8)
    thr2_p,  = ax2.plot([], [], color="tab:gray", linestyle=":", linewidth=0.9, alpha=0.7, label="+/-thr")
    thr2_n,  = ax2.plot([], [], color="tab:gray", linestyle=":", linewidth=0.9, alpha=0.7)
    ax2.axhline(0, color="black", linewidth=0.5, alpha=0.4)
    ax2.set_ylabel("S2 bandpass (mV, raw)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right", fontsize=8)

    # Panel 3: verified-breath timeline. Each green line marks a
    # breath verified by BOTH sensors; each orange line marks a fast-
    # mode breath verified by S1 alone. Shaded regions bracket each
    # breath CYCLE - the span from one verified breath to the next.
    # The two normalised waveforms are drawn faintly underneath as
    # context (showing the phase shift between strain and thermal),
    # but the breath logic uses the per-sensor events, not these
    # traces.
    line_n1, = ax3.plot([], [], color="tab:green",  linewidth=0.7, alpha=0.30, label="S1 norm.")
    line_n2, = ax3.plot([], [], color="tab:purple", linewidth=0.7, alpha=0.30, label="S2 norm.")
    ax3.axhline(0, color="black", linewidth=0.5, alpha=0.4)
    ax3.set_ylabel("Verified breaths")
    ax3.set_xlabel("Time (s)")
    ax3.set_title("Verified breaths  |  green = S2 sourced (slow / moderate)  |  orange = S1 sourced (fast)  |  shaded = one breath cycle")
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc="upper right", fontsize=8)

    # Per-panel breath markers (S1 and S2 panels): red dashed verticals
    # at each sensor's own peak time. These mark every per-sensor
    # detection; the fused markers in panel 3 are the merged version.
    breath_lines1 = LineCollection([], colors="tab:red", linestyles="--",
                                   linewidth=1.0, alpha=0.7)
    breath_lines2 = LineCollection([], colors="tab:red", linestyles="--",
                                   linewidth=1.0, alpha=0.7)
    ax1.add_collection(breath_lines1)
    ax2.add_collection(breath_lines2)

    # Fused panel: separate collections for paired vs fast-mode
    # markers, plus the translucent spans BETWEEN consecutive
    # verified breaths (each span = one breath cycle).
    breath_paired   = LineCollection([], colors="tab:green",  linestyles="-",
                                     linewidth=1.6, alpha=0.95)
    breath_fastmode = LineCollection([], colors="tab:orange", linestyles="-",
                                     linewidth=1.6, alpha=0.95)
    cycle_spans     = PolyCollection([], facecolors="tab:green",
                                     edgecolors="none", alpha=0.12)
    ax3.add_collection(cycle_spans)
    ax3.add_collection(breath_paired)
    ax3.add_collection(breath_fastmode)
    plt.tight_layout()

    # Visible-window buffers (raw bandpass in mV; normalised for panel 3).
    ts             = []
    s1_mv, s2_mv   = [], []
    s1_norm_buf, s2_norm_buf = [], []

    v_scale  = VREF_V / ADC_FULLSCALE
    mv_scale = v_scale * 1000.0

    t0 = time.time()
    n  = 0
    last_thr1 = S1_FLOOR
    last_thr2 = S2_FLOOR

    # One-shot flag for the post-warmup beep below.
    warmup_beep_done = False

    print("[+] Streaming S1+S2 fusion. Ctrl+C or close window to stop.\n")
    try:
        while plt.fignum_exists(fig.number):
            line = ser.readline().decode("utf-8", errors="ignore").strip()
            if not line.startswith("DATA"):
                continue
            fields = {kv.split(":")[0]: kv.split(":")[1]
                      for kv in line.split("|")[1:] if ":" in kv}
            if "S1" not in fields or "S2" not in fields:
                continue
            try:
                y1 = int(fields["S1"])
                y2 = int(fields["S2"])
            except ValueError:
                continue
            now_t = time.time() - t0

            # One-shot audible cue at the moment warm-up ends. From this
            # point on the SensorChain detectors are allowed to emit
            # events, so the beep tells the user "the script is now
            # actually listening for breaths." Non-blocking via a daemon
            # thread; falls back to silent on non-Windows.
            if (not warmup_beep_done) and now_t >= WARMUP_S:
                warmup_beep_done = True
                if winsound is not None:
                    import threading as _th
                    _th.Thread(target=winsound.Beep,
                               args=(1000, 150), daemon=True).start()

            # Scale each sensor's adaptive threshold based on the
            # fuser's current mode, so the trusted sensor stays
            # sensitive to the (shallower) rapid breaths in its mode:
            #   slow     : S1 1.0, S2 1.0
            #   moderate : S1 1.0, S2 scaled  (S2 still trusted)
            #   fast     : S1 scaled, S2 1.0  (S1 takes over)
            mode = fuser.current_mode()
            s1.threshold_scale = FAST_MODE_S1_THR_SCALE     if mode == "fast"     else 1.0
            s2.threshold_scale = MODERATE_MODE_S2_THR_SCALE if mode == "moderate" else 1.0

            bp1, last_thr1, evt1 = s1.step(y1, now_t)
            bp2, last_thr2, evt2 = s2.step(y2, now_t)
            # Drive the fusion state machine. Both ingest() (when an
            # S2 corroborates an S1, or when fast mode emits per-S1)
            # and tick() can return a verified-breath event; collect
            # them all and print on the same code path below.
            fused_events = []
            if evt1 is not None:
                e = fuser.ingest("S1", evt1[0], evt1[1])
                if e is not None:
                    fused_events.append(e)
            if evt2 is not None:
                e = fuser.ingest("S2", evt2[0], evt2[1])
                if e is not None:
                    fused_events.append(e)

            # Envelope-normalise each channel for panel 3 only; the
            # phase shift between the two traces is what makes the
            # "we cannot fuse waveforms" point visually.
            norm1 = bp1 / max(s1.envelope, S1_FLOOR)
            norm2 = bp2 / max(s2.envelope, S2_FLOOR)

            ts.append(now_t)
            s1_mv.append(bp1 * mv_scale)
            # Plot S2 in its NATURAL / raw polarity so the visual
            # matches what the bridge actually outputs (inhale =
            # trough, exhale = peak). The detector still works on
            # the inverted bp internally, so peak_time corresponds
            # to the inverted-peak == raw-TROUGH == end of inhale.
            # Drawing the red breath marker at peak_time therefore
            # lands on the trough here, lining up visually with the
            # S1 peak (which is also end of inhale).
            s2_plot_mv = (-bp2 if S2_INVERT else bp2) * mv_scale
            s2_mv.append(s2_plot_mv)
            s1_norm_buf.append(norm1)
            s2_norm_buf.append(norm2)
            n += 1

            # One console line per completed breath CYCLE - matches the
            # shaded green band on the fused panel one-for-one. A cycle
            # spans the previous verified-breath time to the current
            # one, so the very first verified-breath event only opens
            # cycle 1 and prints nothing; the second event closes it
            # and is reported as "cycle #1 start=... end=...". Tag
            # indicates which sensor sourced the breath that CLOSED
            # the cycle: (S2) for slow/moderate, (S1) for fast.
            for fused_evt in fused_events:
                idx, peak_t, _, interval, bpm, contrib = fused_evt
                tag = ",".join(contrib)
                if interval is None:
                    # First verified breath: marks the START of cycle 1;
                    # nothing to report until cycle 1 closes on the next
                    # verified breath.
                    print(f"[cycle #1 open]  start={peak_t:6.2f}s  ({tag})")
                else:
                    cycle_no = idx - 1
                    start_t  = peak_t - interval
                    print(f"[cycle #{cycle_no}]  start={start_t:6.2f}s  "
                          f"end={peak_t:6.2f}s  duration={interval:5.2f}s  "
                          f"-> {bpm:5.1f} BPM  ({tag})")

            # --- Window slide -----------------------------------------
            cutoff = ts[-1] - WINDOW_S
            drop = 0
            for t in ts:
                if t < cutoff: drop += 1
                else: break
            if drop:
                del ts[:drop]; del s1_mv[:drop]; del s2_mv[:drop]
                del s1_norm_buf[:drop]; del s2_norm_buf[:drop]
            while s1.breath_times and s1.breath_times[0] < cutoff:
                s1.breath_times.pop(0)
            while s2.breath_times and s2.breath_times[0] < cutoff:
                s2.breath_times.pop(0)
            while fuser.cycles_start and fuser.cycles_start[0] < cutoff:
                fuser.cycles_start.pop(0)
                fuser.cycles_end.pop(0)
                fuser.cycles_paired.pop(0)

            # --- Redraw -----------------------------------------------
            if n % UPDATE_EVERY == 0:
                line_s1.set_data(ts, s1_mv)
                line_s2.set_data(ts, s2_mv)
                line_n1.set_data(ts, s1_norm_buf)
                line_n2.set_data(ts, s2_norm_buf)
                thr1_mv = last_thr1 * mv_scale
                thr2_mv = last_thr2 * mv_scale
                thr1_p.set_data([ts[0], ts[-1]], [+thr1_mv, +thr1_mv])
                thr1_n.set_data([ts[0], ts[-1]], [-thr1_mv, -thr1_mv])
                thr2_p.set_data([ts[0], ts[-1]], [+thr2_mv, +thr2_mv])
                thr2_n.set_data([ts[0], ts[-1]], [-thr2_mv, -thr2_mv])

                ax1.set_xlim(ts[0], ts[-1] if ts[-1] > ts[0] else ts[0] + 1e-3)
                for ax in (ax1, ax2, ax3):
                    ax.relim()
                    ax.autoscale_view(scalex=False, scaley=True)

                y1l = ax1.get_ylim(); y2l = ax2.get_ylim(); y3l = ax3.get_ylim()
                breath_lines1.set_segments([[(x, y1l[0]), (x, y1l[1])] for x in s1.breath_times])
                breath_lines2.set_segments([[(x, y2l[0]), (x, y2l[1])] for x in s2.breath_times])
                # Split verified breaths into S1+S2 (green) and fast-
                # mode S1-only (orange) marker sets.
                paired_x   = [t for t, p in zip(fuser.cycles_start, fuser.cycles_paired) if p]
                fastmode_x = [t for t, p in zip(fuser.cycles_start, fuser.cycles_paired) if not p]
                breath_paired.set_segments([[(x, y3l[0]), (x, y3l[1])] for x in paired_x])
                breath_fastmode.set_segments([[(x, y3l[0]), (x, y3l[1])] for x in fastmode_x])
                # Span from each verified breath to the next = one
                # breath cycle. (Last verified breath has no
                # following cycle yet, so no trailing span.)
                spans = []
                for i in range(len(fuser.cycles_start) - 1):
                    st_t = fuser.cycles_start[i]
                    en_t = fuser.cycles_start[i + 1]
                    spans.append([(st_t, y3l[0]), (st_t, y3l[1]),
                                  (en_t, y3l[1]), (en_t, y3l[0])])
                cycle_spans.set_verts(spans)

                fig.canvas.draw_idle()
                fig.canvas.flush_events()
    except KeyboardInterrupt:
        print("\n[+] Stopped by user.")
    finally:
        ser.close()
        plt.ioff()
        if plt.fignum_exists(fig.number):
            plt.show()


if __name__ == "__main__":
    main()
