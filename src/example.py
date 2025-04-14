from spectool.spectool import load_spectrum, determine_continuum, normalize_spectrum, determine_fwhm
file = '/home/jbento/code/ljmu-fwhm-tool/data/arcturus.txt'

wavelength, flux = load_spectrum(file)

continuum = determine_continuum(flux)

normalized_flux = normalize_spectrum(flux, continuum)
results = determine_fwhm(wavelength, normalized_flux)