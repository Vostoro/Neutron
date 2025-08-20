import numpy as np
from math import log, sqrt
from nuclei_data import data, Na, r_0, m0c2, e, me
import photoPB_data
from photoPB_data import photo_pb_E


class Nucleus(object):
    def __init__(self, Ne, E0=0,  element="Fe", uplink=True, is_electron=True):
        self.E1 = 0.0253
        self.E0 = E0
        self.energy_widening()
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
            self.uplinking_cs()
    def energy_widening(self):
        resulting_electron_spectrum = []
        self.electron_energy_indexes = []
        for i in range(0, len(self.FULL_ENERGY_GRID)):
            power = -((self.FULL_ENERGY_GRID[i] - self.E0) / (self.E0/(10**3))) ** 2
            temp = self.Ne * pow(2.71828, power)
            if temp >= 0:
                resulting_electron_spectrum.append(temp)
                self.electron_energy_indexes.append(i)
        resulting_electron_spectrum = np.array(resulting_electron_spectrum)
        self.Ne = resulting_electron_spectrum



    def recalculate_crosssections(self, steps=pow(10, 5), filename="interpvalues.txt"):
        print('Recalcultaing')
        new_file = open(filename, "w")
        N = steps  # Number of steps
        Energies = np.zeros(N)  # Photon energy grid
        Brems_CS = np.zeros(N)  # Bremsstrahlung CS
        E0 = 100 * pow(10, 6)  # Max energy in eV
        E1 = self.E1 / m0c2  # Changing both energies to m0c2 units
        E0 = E0 / m0c2

        for i in range(1, N):
            Energies[i] = E1 + i * (E0 - E1) / N
            E = E0 - Energies[i]
            Brems_CS[i] = 4 * pow(self.Z, 2) * pow(r_0, 2) / 137 * Energies[i] * (
                    (1 + pow(E / E0, 2) - 2 / 3 * E / E0) * log(183 * pow(self.Z, -1)) + 1 / 9 * E / E0)
        # Changing energy back to eV and saving grid to file

        new_file.write("Energy grid\n")
        for i in range(0, len(Energies)):
            Energies[i] = Energies[i] * m0c2
            new_file.write(str(Energies[i]) + ', ')
        new_file.write('\n')

        new_file.write("Bremsstrahlung cs grid\n")
        for sigma in Brems_CS:
            new_file.write(str(sigma) + ', ')
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
            Photoshielding_CS[i] = np.interp(Energies[i], photoPB_data.photo_pb_E, photoPB_data.photo_pb_cs, left=0,
                                             right=0)
        new_file.write("PB Slowing cs grid\n")
        for cross in Photoshielding_CS:
            new_file.write(str(cross) + ', ')
        new_file.write("\n")

        NEUTRON_ENERGY_GRID = np.zeros(N)
        for i in range(0, len(NEUTRON_ENERGY_GRID)):
            if (Energies[i] - 11.1979) > 0:
                NEUTRON_ENERGY_GRID[i] = 55 / 56 * (Energies[i] - 11.1979)
        '''
                NEUTRON_ENERGY_GRID = np.zeros(pow(10,6))
                for i in range(0,len(NEUTRON_ENERGY_GRID)):
                    Eph = Energies[i]
                    mn = 939.565
                    a = 56/55
                    b = 2*Eph/(55*mn)
                    c = (Eph**2)/(55*(mn**2)) - 2*Eph/mn
                    D = -(b**2 - 4*a*c)
                    #print((-b + sqrt(D))/(2*a))
                    if D > 0:
                        NEUTRON_ENERGY_GRID[i] = 1000*(-b + sqrt(D))/(2*a)

        '''

        # Neutron multiplication cross-section
        Neutron_n2n_CS = np.zeros(N)
        Neutron_n3n_CS = np.zeros(N)
        for i in range(1, N):
            Neutron_n2n_CS[i] = np.interp(NEUTRON_ENERGY_GRID[i], photoPB_data.n_pb_2n_E, photoPB_data.n_pb_2n_cs,
                                          left=0, right=0)
            Neutron_n3n_CS[i] = np.interp(NEUTRON_ENERGY_GRID[i], photoPB_data.n_pb_3n_E, photoPB_data.n_pb_3n_cs,
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
        Neutron_detection_cs = np.interp(0.00253, photoPB_data.n_B10_alpha_E, photoPB_data.n_B10_alpha_cs, left=0,
                                         right=0)
        new_file.write("thermal neutron detection cross section\n")
        new_file.write(str(Neutron_detection_cs) + "\n")

        # Photon detection cross-section
        Photon_detection_cs = np.zeros(N)
        for i in range(1, N):
            Photon_detection_cs[i] = np.interp(Energies[i], self.He_E, self.He_cs, left=0, right=0)

        new_file.write("Photon detection cross section\n")
        for cross in Photon_detection_cs:
            new_file.write(str(cross) + ', ')
        new_file.write("\n")
        new_file.close()

    def uplinking_cs(self):
        opened = open("interpvalues.txt")
        file = opened.readlines()
        N = 10 ** 5

        self.FULL_ENERGY_GRID = np.zeros(N)
        string = (file[1].strip('\n')).split(', ')
        for i in range(0, len(self.FULL_ENERGY_GRID)):
            if string[i] != '':
                self.FULL_ENERGY_GRID[i] = float(string[i])

        self.BREMSSTRAHLUNG_GRID = np.zeros(N)
        string = (file[3].strip('\n')).split(', ')
        Total_brems = 0
        for i in range(0, len(self.BREMSSTRAHLUNG_GRID)):
            if string[i] != '':
                self.BREMSSTRAHLUNG_GRID[N - i - 1] = float(string[i]) * 10 ** -12
                Total_brems += float(string[i])



        self.BREMSSTRAHLUNG_GRID = Brems_CS

        self.PHOTONEUTRON_GRID = np.zeros(N)
        string = (file[5].strip('\n')).split(', ')
        for i in range(0, len(self.PHOTONEUTRON_GRID)):
            if string[i] != '':
                self.PHOTONEUTRON_GRID[i] = float(string[i])

        self.PHOTON_SLOWING_GRID = np.zeros(N)
        string = (file[7].strip('\n')).split(', ')
        for i in range(0, len(self.PHOTON_SLOWING_GRID)):
            if string[i] != '':
                self.PHOTON_SLOWING_GRID[i] = float(string[i])

        NEUTRON_ENERGY_GRID = np.zeros(N)
        string = (file[9].strip('\n')).split(', ')
        for i in range(0, len(NEUTRON_ENERGY_GRID)):
            if string[i] != '':
                NEUTRON_ENERGY_GRID[i] = float(string[i])

        self.NEUTRON_N2N_GRID = np.zeros(N)
        string = (file[11].strip('\n')).split(', ')
        for i in range(0, len(self.NEUTRON_N2N_GRID)):
            if string[i] != '':
                self.NEUTRON_N2N_GRID[i] = float(string[i])
        # if self.NEUTRON_N2N_GRID[i] != 0:
        #	print(f'Not a zero, i = {i}')

        self.NEUTRON_N3N_GRID = np.zeros(N)
        string = (file[13].strip('\n')).split(', ')
        for i in range(0, len(self.NEUTRON_N3N_GRID)):
            if string[i] != '':
                self.NEUTRON_N3N_GRID[i] = float(string[i])
        # if self.NEUTRON_N3N_GRID[i] != 0:
        #	print(f'Not a zero, i = {i}')

        string = (file[15].strip('\n'))
        NEUTRON_DETECTION = float(string)

        self.PHOTON_DETECTION_GRID = np.zeros(N)
        string = (file[17].strip('\n')).split(', ')
        if not self.is_electron:
            for i in range(0, len(self.PHOTON_DETECTION_GRID)):
                if string[i] != '':
                    self.PHOTON_DETECTION_GRID[i] = float(string[i]) * pow(10, -8)
        else:
            for i in range(0, len(self.PHOTON_DETECTION_GRID)):
                if string[i] != '':
                    self.PHOTON_DETECTION_GRID[i] = float(string[i])
        Flag = False
        i = 0
        # Using only interpolated values up to the max electron energy
        self.PhotonEnergies = []
        while not Flag:
            try:
                if type(self.E0).__module__ == np.__name__:
                    if self.FULL_ENERGY_GRID[i] <= np.max(self.E0):
                        self.PhotonEnergies.append(self.FULL_ENERGY_GRID[i])
                    else:
                        Flag = True
                else:
                    if self.FULL_ENERGY_GRID[i] <= self.E0:
                        self.PhotonEnergies.append(self.FULL_ENERGY_GRID[i])
                    else:
                        Flag = True
                i += 1
            except IndexError:
                break

        self.Neutron_detection = NEUTRON_DETECTION
        if self.is_electron:
            self.PhotonEnergies = np.array(self.PhotonEnergies)
            # Saving length of the energy array to save time on using the 'len' function
            self.PhotonEnergies_length = len(self.PhotonEnergies)
            self.Bremsstrahlung_cs = self.BREMSSTRAHLUNG_GRID[:self.PhotonEnergies_length]
            self.Photoneutron_cs = self.PHOTONEUTRON_GRID[:self.PhotonEnergies_length]
            self.PhotonSlowing_cs = self.PHOTON_SLOWING_GRID[:self.PhotonEnergies_length]
            self.Neutron_Energies = NEUTRON_ENERGY_GRID[:self.PhotonEnergies_length]
            self.n2n_cs = self.NEUTRON_N2N_GRID[:self.PhotonEnergies_length]
            self.n3n_cs = self.NEUTRON_N3N_GRID[:self.PhotonEnergies_length]
            self.Photon_detection = self.PHOTON_DETECTION_GRID[:self.PhotonEnergies_length]
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
        Nph = np.zeros(self.PhotonEnergies_length)
        Energies = np.zeros(N)  # Photon energy grid
        Brems_CS = np.zeros(N)  # Bremsstrahlung CS
        E0 = self.E0  # Electron energy in eV
        E1 = self.E1 / m0c2  # Changing both energies to m0c2 units
        E0 = E0 / m0c2

        for j in range(1, len(self.Ne)):
            Brems_CS = np.zeros(self.PhotonEnergies_length)
            E0 = self.FULL_ENERGY_GRID[self.electron_energy_indexes[j]]  # Electron energy in eV
            E1 = self.E1 / m0c2  # Changing both energies to m0c2 units
            E0 = E0 / m0c2
            Energies[i] = E1 + i * (E0 - E1) / N
            E = E0 - Energies[i]
            Brems_CS[i] = 4 * pow(self.Z, 2) * pow(r_0, 2) / 137 * Energies[i] * (
                    (1 + pow(E / E0, 2) - 2 / 3 * E / E0) * log(183 * pow(self.Z, -1)) + 1 / 9 * E / E0)
            if filtering == 0:
                for i in range(1, self.PhotonEnergies_length):
                    Nph[i] = self.Bremsstrahlung_cs[i] * self.rho * self.Ne
            else:
                for j in range(0, len(self.Ne)):
                    for i in range(1, self.PhotonEnergies_length):
                        koef = 1 - (2.71828 ** (- filtering * self.PhotonEnergies[i]))
                        Nph[i] = self.Ne[j] * self.Bremsstrahlung_cs[i] * koef * self.rho

        self.BREMSSTRAHLUNG_GRID = None
        self.Bremsstrahlung_cs = None
        return Nph

    def calculate_neutrons(self, Nph, l_wall, indices=[]):
        total = 0
        Nn = np.zeros(len(Nph))
        if self.is_electron:
            for i in range(1, len(Nph)):
                Nn[i] = Nph[i] * (1 - 2.71828 ** (- l_wall * self.Photoneutron_cs[i] * self.rho))
            self.PHOTONEUTRON_GRID = None
            self.Photoneutron_cs = None
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