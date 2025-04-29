###################################
# bespoke bandpass filter in base R 
###################################

# The following function generates the frequency axis corresponding to the output of a Fast Fourier Transform (FFT)
fftfreq <- function(n, d) {
  # Generate frequency axis corresponding to FFT output
  if (n %% 2 == 0) {
    # For even number of data points
    return(seq(0, 1, length.out = n/2 + 1) / d)
  } else {
    # For odd number of data points
    return(seq(0, 1, length.out = (n-1)/2) / d)
  }
}

###################################
bpass <- function(signal, time, period_low, period_high, padfac = 2, win = "none", alpha = 3) {
  # Compute the length of the signal
  n <- length(signal)
  # Calculate sampling interval from time vector
  d <- time[2] - time[1]
  # Remove the mean
  signal <- signal - mean(signal)
  # Zero-padding
  pad_length <- n * padfac
  padded_signal <- c(signal, rep(0, pad_length - n))
  # Compute the FFT of the padded signal
  signal_fft <- fft(padded_signal)
  # Generate the frequency axis using the fftfreq function above
  freq_axis <- fftfreq(length(padded_signal), d)
  # Convert periods to frequencies
  freq_low <- 1 / period_low
  freq_high <- 1 / period_high
  # Create a logical vector indicating frequencies within the bandpass range
  in_bandpass <- (abs(freq_axis) >= freq_low) & (abs(freq_axis) <= freq_high)
  # Apply tapering if specified
  if (win == "gaussian") {
    center_freq <- (freq_low + freq_high) / 2
    bw <- (freq_high - freq_low) / 2
    taper <- exp(-0.5 * ((abs(freq_axis) - center_freq) / bw) ^ 2 * alpha)
    signal_fft <- signal_fft * taper
  } else if (win == "cosine") {
    taper <- rep(1, length(signal_fft))
    p <- 0.25
    bpts <- sum(in_bandpass)
    Pbpts <- round(bpts * p * 0.5)
    for (i in 1:Pbpts) {
      taper[i] <- 0.5 * (1 - cos((2 * pi * i) / (2 * Pbpts + 1)))
      taper[bpts + 1 - i] <- taper[i]
    }
    signal_fft[in_bandpass] <- signal_fft[in_bandpass] * taper
  }
  # Zero out frequencies outside the bandpass range
  signal_fft[!in_bandpass] <- 0
  # Compute the inverse FFT to get the filtered signal
  filtered_signal <- Re(fft(signal_fft, inverse = TRUE) / length(padded_signal))
  # Remove the padding and add the mean back
  filtered_signal <- filtered_signal[1:n] + mean(signal)
  return(filtered_signal)
}
###################################
