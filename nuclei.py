import numpy as np
from math import log, sqrt
from nuclei_data import *

class Nucleus(object):
    def __init__(self, Ne, E0,  element="Fe", uplink=True, is_electron=True):
        self.E1 = 0.0253
        self.E0 = E0
        self.element = element
        self.M = data[self.element]["M"]
        self.Z = data[self.element]["Z"]
        self.rho = data[self.element]["rho"] * 1000000 * Na / data[self.element]["M"]
        self.filename = data[self.element]["filename"]
        self.energymod = data[self.element]["energymod"]
        self.crossectionmod = data[self.element]["crosssectionmod"]
        self.Pb_rho = 11.35 * 1000000 * Na / 208
        self.is_electron = is_electron
        R = 8.31
        T = 300
        p = 1 * 10 ^ 5
        He_rho = (p * 3) / (R * T)
        self.Ne = Ne
        # He_rho = 0.1346
        self.He_rho = He_rho * 1000000 * Na / 3

        self.B10_rho = 2.34 * 1000000 * Na / 10

        opened = open(self.filename)
        file = opened.readlines()
        E = np.zeros([len(file)])
        cs = np.zeros([len(file)])
        i = 0
        for line in file:
            line = line.replace("\n", "")
            E[i] = float(line.split()[0]) * self.energymod
            cs[i] = float(line.split()[1]) * self.crossectionmod
            i += 1
        opened.close()

        self.E = E
        self.cs = cs

        opened = open("data/Helium_Photon_detection_Sig.txt")
        file = opened.readlines()
        self.He_E = np.zeros([len(file)])
        self.He_cs = np.zeros([len(file)])
        i = 0
        for line in file:
            line = line.replace("\n", "")
            self.He_E[i] = float(line.split()[0]) * pow(10, 6)
            self.He_cs[i] = (float(line.split()[1]) + float(line.split()[2]) + float(line.split()[3])) * pow(10, -28)
            i += 1
        opened.close()

        if uplink:
            self.energy_widening()
            self.uplinking_cs()

    def energy_widening(self, sigma_fraction=1e-3, sigma_eV=None, sigma_multiplier=5,
                        threshold_frac=1e-8):
        """
        Create a (normalized) electron-energy distribution on FULL_ENERGY_GRID
        representing a Gaussian spread around self.E0.

        Parameters
        ----------
        sigma_fraction : float
            If sigma_eV is None, sigma = sigma_fraction * E0.
        sigma_eV : float or None
            If provided, use this absolute sigma (in eV) instead of sigma_fraction.
        sigma_multiplier : float
            Keep only grid points within +/- sigma_multiplier * sigma (for speed).
        threshold_frac : float
            Bins with population < threshold_frac * total_electrons will be omitted
            from self.electron_energy_indexes.

        Effects
        -------
        Sets self.Ne to the resulting array (counts per energy bin) and
        self.electron_energy_indexes to the list of indices with significant counts.

        Returns
        -------
        resulting_electron_spectrum : np.ndarray
            The array of electron counts per FULL_ENERGY_GRID bin.
        """

        opened = open("interpvalues.txt")
        file = opened.readlines()
        self.amount_of_elements = int(file[0].strip('\n'))

        self.FULL_ENERGY_GRID = np.zeros(self.amount_of_elements)
        string = (file[2].strip('\n')).split(', ')
        for i in range(0, len(self.FULL_ENERGY_GRID)):
            if string[i] != '':
                self.FULL_ENERGY_GRID[i] = float(string[i])


        Egrid = np.asarray(self.FULL_ENERGY_GRID, dtype=float)
        if Egrid.size == 0:
            raise ValueError("FULL_ENERGY_GRID is empty")

        E0 = float(self.E0)

        # Determine total number of electrons (if self.Ne is scalar or already an array)
        Ne_attr = np.asarray(self.Ne)
        if Ne_attr.ndim == 0:
            total_electrons = float(Ne_attr)
        else:
            total_electrons = float(np.sum(Ne_attr))

        if sigma_eV is None:
            sigma = max(abs(E0) * float(sigma_fraction), 1e-12)
        else:
            sigma = float(sigma_eV)
            if sigma <= 0:
                raise ValueError("sigma_eV must be positive")

        # Build Gaussian (standard form: exp(-0.5 * ((x/sigma)^2)))
        x = (Egrid - E0) / sigma
        gauss = np.exp(-0.5 * x ** 2)

        # Optionally zero out far tails for efficiency
        mask_range = np.abs(x) <= sigma_multiplier
        gauss = gauss * mask_range

        # If everything masked out (e.g., sigma very small), place all electrons in nearest bin
        if gauss.sum() <= 0:
            spectrum = np.zeros_like(gauss)
            idx = np.abs(Egrid - E0).argmin()
            spectrum[idx] = total_electrons
        else:
            spectrum = gauss / gauss.sum() * total_electrons

        # Build electron_energy_indexes: only bins above threshold
        threshold = threshold_frac * total_electrons
        idxs = np.nonzero(spectrum > threshold)[0].tolist()

        # Save to object (overwrites self.Ne as an array of bin counts)
        self.electron_energy_indexes = idxs
        self.Ne = spectrum

        return spectrum

    def recalculate_crosssections(self, steps=pow(10, 5), filename="interpvalues.txt"):
        """
        Recalculate cross-section by interpolating the values from nuclei_data.

        Parameters
        ----------
        steps : int
            The amount of steps between 0.0253 and max energy
        filename : String
            The name of the file in which the interpolated values are saved.

        Effects
        -------
        Creates a file with values of energy grid and cross-sections.

        Returns
        -------
        Nothing

        """

        new_file = open(filename, "w")
        N = steps  # Number of steps
        Energies = np.zeros(N)  # Photon energy grid
        new_file.write(f"{steps}\n")
        dE = (self.E0 - self.E1)/N
        # Saving grid to file
        new_file.write("Energy grid\n")
        for i in range(0, N):
            Energies[i] = self.E1 + i * dE
            new_file.write(str(Energies[i]) + ', ')
        new_file.write('\n')

        new_file.write("LEGACY LINES\n")
        new_file.write("\n")

        # Photoneutron cross-section
        Photoneut_CS = np.zeros(N)
        for i in range(1, N):
            Photoneut_CS[i] = np.interp(Energies[i], self.E, self.cs, left=0, right=0)

        new_file.write("Photoneutron cs grid\n")
        for cross in Photoneut_CS:
            new_file.write(str(cross) + ', ')
        new_file.write("\n")

        # Photon shielding cross-section
        Photoshielding_CS = np.zeros(N)
        for i in range(1, N):
            Photoshielding_CS[i] = np.interp(Energies[i], photo_pb_E, photo_pb_cs, left=0,
                                             right=0)
        new_file.write("PB Slowing cs grid\n")
        for cross in Photoshielding_CS:
            new_file.write(str(cross) + ', ')
        new_file.write("\n")

        NEUTRON_ENERGY_GRID = np.zeros(N)
        for i in range(0, len(NEUTRON_ENERGY_GRID)):
            if (Energies[i] - 11.1979) > 0:
                NEUTRON_ENERGY_GRID[i] = 55 / 56 * (Energies[i] - 11.1979)

        # Neutron multiplication cross-section
        Neutron_n2n_CS = np.zeros(N)
        Neutron_n3n_CS = np.zeros(N)
        for i in range(1, N):
            Neutron_n2n_CS[i] = np.interp(NEUTRON_ENERGY_GRID[i], n_pb_2n_E, n_pb_2n_cs,
                                          left=0, right=0)
            Neutron_n3n_CS[i] = np.interp(NEUTRON_ENERGY_GRID[i], n_pb_3n_E, n_pb_3n_cs,
                                          left=0, right=0)

        # Neutron spectrum differs from the photon spectrum and thus its energies are written separately
        new_file.write("photoneutron energy spectrum grid\n")
        for energy in NEUTRON_ENERGY_GRID:
            new_file.write(str(energy) + ', ')
        new_file.write("\n")
        new_file.write("n2n cs grid\n")
        for cross in Neutron_n2n_CS:
            new_file.write(str(cross) + ', ')
        new_file.write("\n")
        new_file.write("n3n cs grid\n")
        for cross in Neutron_n3n_CS:
            new_file.write(str(cross) + ', ')
        new_file.write("\n")

        # Neutron detection cross-section
        Neutron_detection_cs = np.interp(0.00253, n_B10_alpha_E, n_B10_alpha_cs, left=0,
                                         right=0)
        new_file.write("thermal neutron detection cross section grid\n")
        new_file.write(str(Neutron_detection_cs) + "\n")

        # Photon detection cross-section
        Photon_detection_cs = np.zeros(N)
        for i in range(1, N):
            Photon_detection_cs[i] = np.interp(Energies[i], self.He_E, self.He_cs, left=0, right=0)

        new_file.write("Photon detection cross section grid\n")
        for cross in Photon_detection_cs:
            new_file.write(str(cross) + ', ')
        new_file.write("\n")
        new_file.close()

    def uplinking_cs(self, filename="interpvalues.txt"):
        """
                Gets and sets up the energy grid and cross-section grid from file and chooses only the energies lower
                than the starting energy of electrons.

                Parameters
                ----------
                filename : String
                    The name of the file from which the interpolated values read.

                Effects
                -------
                Sets up energy grid and cross-cross sections. Saves only values up to  max value of self.E0
                    self.PHOTONEUTRON_GRID - photoneutron cross-section generation for material chosen at the last recalculation
                    self.PHOTON_SLOWING_GRID - photon slowing cross-section for Pb
                    self.NEUTRON_N2N_GRID - (n,2n) multiplication cross-section for Pb
                    self.NEUTRON_N3N_GRID - (n,3n) multiplication cross-section for Pb
                    self.PHOTON_DETECTION_GRID - electron pair production cross-section for Helium
                    self.PhotonEnergies - the energy grid
                    self.PhotonEnergies_length - the length of elements in the grid

                Returns
                -------
                Nothing

                """

        opened = open(filename)
        file = opened.readlines()

        self.amount_of_elements = int(file[0].strip('\n'))

        N = self.amount_of_elements

        self.FULL_ENERGY_GRID = np.zeros(N)
        string = (file[2].strip('\n')).split(', ')
        for i in range(0, len(self.FULL_ENERGY_GRID)):
            if string[i] != '':
                self.FULL_ENERGY_GRID[i] = float(string[i])


        PHOTONEUTRON_GRID = np.zeros(N)
        string = (file[6].strip('\n')).split(', ')
        for i in range(0, len(PHOTONEUTRON_GRID)):
            if string[i] != '':
                PHOTONEUTRON_GRID[i] = float(string[i])

        PHOTON_SLOWING_GRID = np.zeros(N)
        string = (file[8].strip('\n')).split(', ')
        for i in range(0, len(PHOTON_SLOWING_GRID)):
            if string[i] != '':
                PHOTON_SLOWING_GRID[i] = float(string[i])

        NEUTRON_ENERGY_GRID = np.zeros(N)
        string = (file[10].strip('\n')).split(', ')
        for i in range(0, len(NEUTRON_ENERGY_GRID)):
            if string[i] != '':
                NEUTRON_ENERGY_GRID[i] = float(string[i])

        NEUTRON_N2N_GRID = np.zeros(N)
        string = (file[12].strip('\n')).split(', ')
        for i in range(0, len(NEUTRON_N2N_GRID)):
            if string[i] != '':
                NEUTRON_N2N_GRID[i] = float(string[i])
        # if self.NEUTRON_N2N_GRID[i] != 0:
        #	print(f'Not a zero, i = {i}')

        NEUTRON_N3N_GRID = np.zeros(N)
        string = (file[14].strip('\n')).split(', ')
        for i in range(0, len(NEUTRON_N3N_GRID)):
            if string[i] != '':
                NEUTRON_N3N_GRID[i] = float(string[i])
        # if self.NEUTRON_N3N_GRID[i] != 0:
        #	print(f'Not a zero, i = {i}')

        string = (file[16].strip('\n'))
        NEUTRON_DETECTION = float(string)

        PHOTON_DETECTION_GRID = np.zeros(N)
        string = (file[18].strip('\n')).split(', ')
        if not self.is_electron:
            for i in range(0, len(PHOTON_DETECTION_GRID)):
                if string[i] != '':
                    PHOTON_DETECTION_GRID[i] = float(string[i]) * pow(10, -8)
        else:
            for i in range(0, len(PHOTON_DETECTION_GRID)):
                if string[i] != '':
                    PHOTON_DETECTION_GRID[i] = float(string[i])


        # Using only interpolated values up to the max electron energy
        Egrid = np.asarray(self.FULL_ENERGY_GRID, dtype=float)
        # cutoff = max(self.E0) works whether self.E0 is scalar or array-like
        cutoff = float(np.max(self.E0))
        tol = 1e-12
        # find rightmost index where Egrid[idx] <= cutoff (use tol to avoid fp issues)
        end_index = np.searchsorted(Egrid, cutoff + tol, side='right')
        # slice from 0 up to end_index (may be empty)
        selected = Egrid[:end_index]
        # store as numpy array (fast for later numeric work)

        self.Neutron_detection = NEUTRON_DETECTION
        if self.is_electron:
            # Saving length of the energy array to save time on using the 'len' function
            self.PhotonEnergies = selected
            self.PhotonEnergies_length = selected.size
            self.Photoneutron_cs = PHOTONEUTRON_GRID[:self.PhotonEnergies_length]
            self.PhotonSlowing_cs = PHOTON_SLOWING_GRID[:self.PhotonEnergies_length]
            self.Neutron_Energies = NEUTRON_ENERGY_GRID[:self.PhotonEnergies_length]
            self.n2n_cs = NEUTRON_N2N_GRID[:self.PhotonEnergies_length]
            self.n3n_cs = NEUTRON_N3N_GRID[:self.PhotonEnergies_length]
            self.Photon_detection = PHOTON_DETECTION_GRID[:self.PhotonEnergies_length]

        else:
            return

    def answer_combing(self, value):
        if type(value) != np.ndarray:
            answer_str = str(round(value))
            crutch = []
            for i in answer_str:
                crutch.append(i)
            answer_str = str(round(int(answer_str) / pow(10, len(crutch) - 1), 2)) + 'e' + str(len(crutch) - 1)
        else:
            answer_str = str(round(sum(value)))
            crutch = []
            for i in answer_str:
                crutch.append(i)
            answer_str = str(round(int(answer_str) / pow(10, len(crutch) - 1), 2)) + 'e' + str(len(crutch) - 1)
        return answer_str

    def preset_photons_spectrum(self, spectrum_type, k, b, c, d):
        Nph = []
        local_energy_indices = []
        if spectrum_type == 'Линейный':
            for i in range(0, len(self.FULL_ENERGY_GRID)):
                temp = self.FULL_ENERGY_GRID[i] * k + b
                if temp >= 0:
                    Nph.append(temp)
                    local_energy_indices.append(i)
        elif spectrum_type == 'Параболический':
            for i in range(0, len(self.FULL_ENERGY_GRID)):
                temp = k * (self.FULL_ENERGY_GRID[i] ** 2) + b * self.FULL_ENERGY_GRID[i] + c
                if temp >= 0:
                    Nph.append(temp)
                    local_energy_indices.append(i)
        elif spectrum_type == 'Кубическая парабола':
            for i in range(0, len(self.FULL_ENERGY_GRID)):
                temp = k * (self.FULL_ENERGY_GRID[i] ** 3) + b * (self.FULL_ENERGY_GRID[i] ** 2) + c * \
                       self.FULL_ENERGY_GRID[i] + d
                if temp >= 0:
                    Nph.append(temp)
                    local_energy_indices.append(i)
        elif spectrum_type == 'Гаусс':
            for i in range(0, len(self.FULL_ENERGY_GRID)):
                power = -((self.FULL_ENERGY_GRID[i] - b) / c) ** 2
                temp = k * pow(2.71828, power) + d
                # temp
                # print(power)
                if temp >= 0:
                    Nph.append(temp)
                    local_energy_indices.append(i)
        elif spectrum_type == 'TEST':
            Nph = 10000
            local_energy_indices = 200000
        Nph = np.array(Nph)
        local_energy_indices = np.array(local_energy_indices)
        return Nph, local_energy_indices

    def calculate_photons(self, filtering=0):
        """
        Calculate number of photons in each photon-energy bin produced by
        bremsstrahlung from electrons in the states described by
        self.electron_energy_indexes / self.FULL_ENERGY_GRID and populations self.Ne.

        Returns
        -------
        Nph : numpy.ndarray
            Photon counts per photon-energy bin (length self.PhotonEnergies_length).
        """
        # Photon energy array in eV
        photon_energies = np.asarray(self.PhotonEnergies)
        n_ph_bins = self.PhotonEnergies_length
        if photon_energies.size != n_ph_bins:
            photon_energies = photon_energies[:n_ph_bins]

        # Prepare output: photon counts per bin
        Nph = np.zeros(n_ph_bins, dtype=float)

        # Convert E1 and any electron energies to m0c2 units for cross-section formula,
        # but we will keep photon_energies in eV for threshold checks and convert when needed.
        E1_eV = float(self.E1)  # in eV
        E1 = E1_eV / m0c2  # dimensionless in units of m0c2

        # Electrons population: allow self.Ne to be scalar or an array
        Ne_arr = np.asarray(self.Ne)
        scalar_Ne = False
        if Ne_arr.ndim == 0:
            scalar_Ne = True

        # Get list of electron-energy indices we'll iterate over.
        # If self.electron_energy_indexes is not present, fall back to a single energy self.E0.
        try:
            electron_indexes = list(self.electron_energy_indexes)
        except Exception:
            # fallback: single energy in self.E0
            electron_indexes = [None]

        # Will store per-electron-state bremsstrahlung cross sections (for debugging / storing)
        brems_cs_grid = []

        # Pre-calc constants appearing in the formula
        Z = float(self.Z)
        rho = float(self.rho)
        prefactor = 4.0 * (Z ** 2) * (r_0 ** 2) / 137.0

        # Loop over electron-energy states
        for idx_i, e_idx in enumerate(electron_indexes):
            # Electron initial energy in eV (before emitting a photon)
            if e_idx is None:
                E0_eV = float(self.E0)
            else:
                E0_eV = float(self.FULL_ENERGY_GRID[e_idx])

            # Convert to m0c2 units for the formula
            E0 = E0_eV / m0c2

            # Maximum photon energy allowed for this electron state:
            # photon_energy <= E0_eV - E1_eV (so final electron energy >= E1_eV)
            max_photon_energy = max(0.0, E0_eV - E1_eV)

            # mask of photon bins that are physically allowed for this electron state
            allowed_mask = photon_energies <= max_photon_energy

            if not np.any(allowed_mask):
                # No allowed photons from this electron energy -> zero contribution
                brems_cs_grid.append(np.zeros(n_ph_bins))
                continue

            # Convert the photon energies to m0c2 units for formula
            photon_m = photon_energies / m0c2

            # Electron energy after emission (E_after = E_before - photon)
            E_after_m = (E0_eV - photon_energies) / m0c2

            # Avoid negative or zero division; keep computation only on allowed_mask
            ratio = np.zeros_like(photon_m)
            ratio[allowed_mask] = E_after_m[allowed_mask] / E0  # E_after / E_before

            # Compute the differential cross-section (vectorized)
            # cs = prefactor * photon_m * ( (1 + ratio^2 - (2/3)*ratio )*log(183/Z) + (1/9)*ratio )
            # Use np.errstate to suppress warnings where ratio is undefined
            with np.errstate(divide='ignore', invalid='ignore'):
                log_term = np.log(183.0 * (Z ** -1.0))
                cs = np.zeros(n_ph_bins, dtype=float)
                cs_expr = (1.0 + ratio ** 2 - (2.0 / 3.0) * ratio) * log_term + (1.0 / 9.0) * ratio
                cs[allowed_mask] = prefactor * photon_m[allowed_mask] * cs_expr[allowed_mask]

            # store per-state cs row
            brems_cs_grid.append(cs)

            # number of electrons in this state
            if scalar_Ne:
                ne_count = float(self.Ne)
            else:
                # assume self.Ne aligns with electron_indexes; if not, we take nearest or 0
                try:
                    ne_count = float(Ne_arr[idx_i])
                except Exception:
                    # fallback: if lengths mismatch, assume zero
                    ne_count = 0.0

            # filtering weight (if filtering==0 => weight 1)
            if filtering == 0:
                weight = 1.0
            else:
                # use exponential attenuation: 1 - exp(- filtering * photon_energy)
                # photon_energies in eV are fine as argument to the user-given filtering factor
                weight = 1.0 - np.exp(- filtering * photon_energies)

            # Add contribution to total photon counts
            # Nph += ne_count * cs * weight * rho
            Nph += ne_count * cs * weight * rho

        return Nph

    def calculate_neutrons(self, Nph, l_wall, indices=[]):
        total = 0
        Nn = np.zeros(len(Nph))
        if self.is_electron:
            for i in range(1, len(Nph)):
                Nn[i] = Nph[i] * (1 - 2.71828 ** (- l_wall * self.Photoneutron_cs[i] * self.rho))
        else:
            for i in range(1, len(Nph)):
                Nn[i] = Nph[i] * (1 - 2.71828 ** (- l_wall * self.PHOTONEUTRON_GRID[indices[i]] * self.rho))
        return Nn

    def calculate_shielded_photons(self, Nph, l_shield, indices=[]):
        output_photons = np.zeros(len(Nph))
        if self.is_electron:
            for i in range(1, len(Nph)):
                output_photons[i] = Nph[i] * (2.71828 ** (- l_shield * self.PhotonSlowing_cs[i] * self.Pb_rho))
            self.PHOTON_SLOWING_GRID = None
            self.PhotonSlowing_cs = None

        else:
            for i in range(1, len(Nph)):
                output_photons[i] = Nph[i] * (
                            2.71828 ** (- l_shield * self.PHOTON_SLOWING_GRID[indices[i]] * self.Pb_rho))
        return output_photons

    def calculate_shielded_neutrons(self, Nn, l_shield, indices=[]):
        output_neutrons = sum(Nn)
        if self.is_electron:
            for i in range(0, len(Nn)):
                n2 = Nn[i] * (1 - 2.71828 ** (- l_shield * self.n2n_cs[i] * self.Pb_rho))
                # print(f'Additionall n2n neutrons = {2*n2}')
                n3 = Nn[i] * (1 - 2.71828 ** (- l_shield * self.n3n_cs[i] * self.Pb_rho))
                output_neutrons += n2 / 2 + 2 * n3
            self.NEUTRON_N2N_GRID = None
            self.NEUTRON_N3N_GRID = None
        else:
            for i in range(0, len(Nn)):
                # print(f'N2n reaction count {(1 - 2.71828 ** (- l_shield * self.NEUTRON_N2N_GRID[indices[i]] * self.Pb_rho))}')
                n2 = Nn[i] * (1 - 2.71828 ** (- l_shield * self.NEUTRON_N2N_GRID[indices[i]] * self.Pb_rho))
                # print(f'Additionall n2n neutrons = {2*n2}')
                n3 = Nn[i] * (1 - 2.71828 ** (- l_shield * self.NEUTRON_N3N_GRID[indices[i]] * self.Pb_rho))
                output_neutrons += n2 / 2 + 2 * n3
        return output_neutrons

    def calculate_signal(self, Nn, Nph, indices=[]):
        thickness = 0.000010
        helium_rad = 0.01
        photons = 0
        output = 0
        if self.is_electron:
            try:
                output = sum(Nn) * (1 - 2.71828 ** (-2 * thickness * self.B10_rho * self.Neutron_detection))
                print(f'inputed neutrons {sum(Nn)}')
                print(f'detected amount = {output / sum(Nn)}')
            except TypeError:
                print(f'inputed neutrons{Nn}')
                output = Nn * (1 - 2.71828 ** (-2 * thickness * self.B10_rho * self.Neutron_detection))
                print(f'detected amount = {output / Nn}')
            print(f'Detected neutron = {output}')
            for i in range(0, len(Nph)):
                if (1 - 2.71828 ** (-helium_rad * self.He_rho * self.Photon_detection[i])) != 0:
                    test = Nph[i] * (1 - 2.71828 ** (-2 * helium_rad * self.He_rho * self.Photon_detection[i]))
                    photons += test
                    output += test
            print(f'Photons = {photons}')
            self.PHOTON_DETECTION_GRID = None
        # self.Photon_detection = None
        else:
            try:
                output = sum(Nn) * (1 - pow(2.71828, (-2 * thickness * self.B10_rho * self.Neutron_detection)))
                print(f'inputed neutrons {sum(Nn)}')
                print(f'detected amount = {output / sum(Nn)}')
            except TypeError:
                print(f'inputed neutrons{Nn}')
                output = Nn * (1 - pow(2.71828, (-2 * thickness * self.B10_rho * self.Neutron_detection)))
                print(f'detected amount = {output / Nn}')
            print(f'Detected neutron = {output}')

            try:
                for i in range(0, len(Nph)):
                    detection = (
                                1 - 2.71828 ** (-2 * helium_rad * self.He_rho * self.PHOTON_DETECTION_GRID[indices[i]]))
                    # print(detection)
                    test = Nph[i] * detection
                    photons += test
                    output += test
                print(f'helium detection {output / sum(Nph)}')
            except TypeError:
                test = Nph * (1 - 2.71828 ** (-2 * helium_rad * self.He_rho * self.PHOTON_DETECTION_GRID[indices]))
                photons = test
                output += test
                print(f'helium detection {output / Nn}')
            print(f'Photons = {photons}')
        return output

    def alloy_calculate_photoneutrons(self, Ne, E0, l=1):
        Fe = Nucleus(element="Fe")
        Fe_n = Fe.calculate_photoneutrons(Ne, E0, l) * 0.72
        Cr = Nucleus(element="Cr")
        Cr_n = Cr.calculate_photoneutrons(Ne, E0, l) * 0.18
        Ni = Nucleus(element="Ni")
        Ni_n = Ni.calculate_photoneutrons(Ne, E0, l) * 0.10
        calculated = Fe_n + Cr_n + Ni_n
        return calculated