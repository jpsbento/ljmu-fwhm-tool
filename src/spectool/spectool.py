
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k
from scipy.optimize import curve_fit

central_wavelengths = {"Halpha": 6563, "Hbeta": 4861, "Hgamma": 4340, "Hdelta": 4102, "Hepsilon": 3970,}

def load_spectrum(file):
    
    wavelength = []
    flux = []
    # The wl of the spectrum is in the first row.
    # The flux is in the second row. 
    # We will read the file line by line and split each line into two values.
    with open(file, 'r') as file:
        for line in file:
            # Process data lines (wavelength and flux values)
            v = line.split()
            if len(v) == 2:  # Ensure there are exactly two values
                wavelength_value = float(v[0])
                flux_value = float(v[1])
                # Filter by max wavelength
                if wavelength_value <= 15000:
                    wavelength.append(wavelength_value)
                    flux.append(flux_value)
    # Convert lists to numpy arrays for easy manipulation
    wavelength = np.array(wavelength)
    flux = np.array(flux)
    return wavelength, flux



def determine_continuum(flux, windowSize=50, threshold=0.8):
    log_flux = np.log(flux)

    continuum = np.zeros_like(flux)


    for i in range(len(flux)):
        half_window = min(i, len(flux) - i - 1, windowSize // 2)
        
        # Define the window range in linear space
        window_start = i - half_window
        window_end = i + half_window + 1
        
        # Calculate the running average in log-flux space, excluding deep absorption lines
        window_flux = log_flux[window_start:window_end]
        filtered_flux = window_flux[window_flux > np.log(threshold * np.median(np.exp(window_flux)))]
        
        # Take the median of the filtered flux, or the median of the whole window if no filtering
        continuum[i] = np.exp(np.median(filtered_flux)) if len(filtered_flux) > 0 else np.exp(np.median(window_flux))

    return continuum

def plot_spectrum(wavelength, flux, continuum):
    plt.figure(figsize=(10, 6))
    plt.plot(wavelength, flux, label="Original Spectrum")
    plt.plot(wavelength, continuum, label=f"COntinuum Fit", linestyle="--")
    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel("Flux")
    plt.legend()
    plt.title("Blackbody Continuum Fit")
    plt.show()



def normalize_spectrum(flux, continuum):
    flux = flux / continuum
    return flux

# Fit Gaussian around each target line
def gaussian(x, amp, mu, sigma):
    return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def determine_fwhm(wavelength, flux, central_wavelengths=central_wavelengths):
    fwhm_results = []


    # Create the inverted spectrum for fitting
    inverted_spectrum = 1 - flux

    plt.title("Spectral Lines with Fitted Gaussians")
    for line, center_wavelength in central_wavelengths.items():
        # Initial guess based on the target line center
        idx = (np.abs(wavelength - center_wavelength)).argmin()
        amp_guess = inverted_spectrum[idx]
        sigma_guess = 10  # Initial guess for sigma
        guess = [amp_guess, center_wavelength, sigma_guess]

        # Narrow region around the line center for fitting
        mask = (wavelength > center_wavelength - 200) & (wavelength < center_wavelength + 200)

        # Check if there are enough points to fit
        if np.sum(mask) < 5:  # Require at least 5 points in the window
            print(f"Not enough data points to fit %s" % line)
        else:   
            popt, _ = curve_fit(gaussian, wavelength[mask], inverted_spectrum[mask], p0=guess)

            # Calculate FWHM from the fit sigma
            amp, mu, sigma = popt
            fwhm = 2 * np.sqrt(2 * np.log(2)) * np.abs(sigma)
            print('The FWHM of %s is %.2f Å' % (line,fwhm))
            fwhm_results.append((line, mu, fwhm, popt, mask))
    
            # Plot the fitted Gaussian, inverted back to negative values
            plt.plot(wavelength[mask], 1-gaussian(wavelength[mask], *popt), linestyle='--', label='%s FWHM=%s Å' % (line,fwhm))
            plt.plot(wavelength[mask], flux[mask], label=line + ' Spectrum', color='black')
    
    plt.show()
    return fwhm_results


def full_analysis(file):
    # Load the spectrum
    wavelength, flux = load_spectrum(file)
    
    # Determine the continuum
    continuum = determine_continuum(flux)
    
    # Normalize the spectrum
    normalized_flux = normalize_spectrum(flux, continuum)
    
    # Plot the original and normalized spectrum
    plot_spectrum(wavelength, flux, continuum)
    
    # Determine FWHM for specified lines
    fwhm_results = determine_fwhm(wavelength, normalized_flux)
    
    return fwhm_results