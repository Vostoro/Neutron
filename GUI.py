import sys, os

import numpy as np
from PyQt6.QtWidgets import QWidget, QMainWindow, QDialog, QDialogButtonBox, QVBoxLayout, QLabel, QApplication, \
    QFileDialog, QStackedWidget, QCheckBox, QComboBox, QPushButton, QTextEdit
from PyQt6 import uic
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

from nuclei import *

energies = {
    "МэВ" : pow(10,6),
    "кэВ": pow(10,3),
    "эВ": 1
}

spectra = {
    "Линейный",
    "Параболический",
    "Гаусс"
}

ui_photon, _ = uic.loadUiType('Elphoneut.ui')

class PhotonWindow(QMainWindow, ui_photon):
    def __init__(self, parent=None):
        super(PhotonWindow, self).__init__(parent)
        self.setupUi(self)
        #self.sensitivity_test()
        #return
        self.press = None
        self.recalculateFlag = False
        self.spectrumPlotFlag = False
        self.patchFlag = False

        self.electron_widgets = [self.amountLabel, self.energyLabel, self.filterCheckBox, self.energyLabel, self.amountLabel, self.amountTextEdit,
                            self.energyTextEdit]
        self.photon_widgets = [self.aLabel, self.bLabel, self.cLabel, self.dLabel, self.aTextEdit,
                               self.bTextEdit, self.cTextEdit, self.dTextEdit,
                               self.spectrumComboBox, self.equationLabel]

        self.methodComboBox.addItem('Расчёт от пучка электронов')
        self.methodComboBox.addItem('Расчёт по спектру фотонов')
        self.methodComboBox.setCurrentText('Расчёт от пучка электронов')
        self.methodChange()

        for spectrum in spectra:
            self.spectrumComboBox.addItem(spectrum)
        for target in data.keys():
            self.targetComboBox.addItem(target)
        self.targetComboBox.setCurrentText("Fe")

        self.recalculateGridCheckstateChanged()

        self.signalFigure = Figure()
        self.signalCanvas = FigureCanvas(self.signalFigure)

        self.spectrumFigure = Figure()
        self.spectrumCanvas = FigureCanvas(self.spectrumFigure)

        self.plot_connect()

        self.recalculateCheckBox.stateChanged.connect(self.recalculateGridCheckstateChanged)
        self.filterCheckBox.stateChanged.connect(self.filterChanged)

        self.methodComboBox.currentTextChanged.connect(self.methodChange)
        self.spectrumComboBox.currentTextChanged.connect(self.spectrumChange)


        self.calculatePushButton.clicked.connect(self.Calculate)

    def beamCalculate(self):
        if self.spectrumPlotFlag:
            self.plotLayout.removeWidget(self.spectrumCanvas)
            self.spectrumPlotFlag = False
        l = self.doubleSpinBox.value() * 100
        l_wall = self.doubleSpinBox_2.value() * 100
        curElement = self.targetComboBox.currentText()
        E0 = self.valueFromText(self.energyTextEdit.toPlainText())
        Ne = self.valueFromText(self.amountTextEdit.toPlainText())
        if self.recalculateFlag:
            curNuclei = Nucleus(E0=self.valueFromText(self.maxTextEdit.toPlainText()), Ne=Ne, element=curElement, uplink=False)
            curNuclei.recalculate_crosssections()
            curNuclei = Nucleus(E0=E0, Ne=Ne, element=curElement, uplink=True, is_electron=True)
        else:
            curNuclei = Nucleus(E0=E0, Ne=Ne, element=curElement, uplink=True, is_electron=True)

        if not self.filterCheckBox.isChecked():
            print(self.valueFromText(self.aTextEdit.toPlainText()))

        electrons = curNuclei.energy_widening()
        electron_energies = curNuclei.FULL_ENERGY_GRID

        plt.plot(electron_energies,electrons)
        plt.xlabel('Энергия электронов, эВ')
        plt.ylabel('Число электронов, шт.')
        plt.show()

        self.photons = curNuclei.calculate_photons(filtering=self.valueFromText(self.aTextEdit.toPlainText()))
        self.spectrumEnergyGrid = curNuclei.PhotonEnergies
        self.plot_spectrum(self.spectrumEnergyGrid, self.photons, is_electron=True)
        neutrons = curNuclei.calculate_neutrons(Nph=self.photons, l_wall=l)
        shielded_photons = curNuclei.calculate_shielded_photons(Nph=self.photons, l_shield=l_wall)
        shielded_neutrons = curNuclei.calculate_shielded_neutrons(Nn=neutrons, l_shield=l_wall)
        self.signal_1 = curNuclei.calculate_signal(Nn=neutrons, Nph=self.photons)
        self.signal_2 = curNuclei.calculate_signal(Nn=shielded_neutrons, Nph=shielded_photons)


    # Calculation from predetermined spectrum
    def spectrumCalculate(self):
        a_coeff = self.valueFromText(self.aTextEdit.toPlainText())
        b_coeff = self.valueFromText(self.bTextEdit.toPlainText())
        c_coeff = self.valueFromText(self.cTextEdit.toPlainText())
        d_coeff = self.valueFromText(self.dTextEdit.toPlainText())
        print(f'read value = {a_coeff}')
        spectrumType = self.spectrumComboBox.currentText()
        l = self.doubleSpinBox.value()
        l_wall = self.doubleSpinBox_2.value()
        curElement = self.targetComboBox.currentText()
        equation = "N = "
        if spectrumType == "Линейный":
            self.equationLabel.setText(f"N = {a_coeff}*E  + {b_coeff}")
        elif spectrumType == "Параболический":
            self.equationLabel.setText(f"N = {a_coeff}*E^2  + {b_coeff}*E + {c_coeff}")
        elif spectrumType == "Кубическая парабола":
            self.equationLabel.setText(f"N = {a_coeff}*E^3  + {b_coeff}*E^2 + {c_coeff}*E + {d_coeff}")
        elif spectrumType == "Гаусс":
            "N = ae^-((E-b)^2/c^2) + d"
            self.equationLabel.setText(f"N = {a_coeff}e^-((-E - {b_coeff})^2/{c_coeff}) + {d_coeff}")

        if self.recalculateFlag:
            print('Starting recalculating')
            curNuclei = Nucleus(element=curElement, uplink=False, is_electron=False)
            curNuclei.recalculate_crosssections()
            print('Starting Uplinking')
            curNuclei = Nucleus(element=curElement, uplink=True, is_electron=False)
            print('Uplink done')
        else:
            curNuclei = Nucleus(E0=0, element=curElement, uplink=True, is_electron=False)
        photons, indexes = curNuclei.preset_photons_spectrum(spectrum_type=spectrumType, k=a_coeff, b=b_coeff,
                                                             c=c_coeff, d=d_coeff)
        neutrons = curNuclei.calculate_neutrons(Nph=photons, l_wall=l, indices=indexes)
        shielded_photons = curNuclei.calculate_shielded_photons(Nph=photons, l_shield=l_wall, indices=indexes)
        shielded_neutrons = curNuclei.calculate_shielded_neutrons(Nn=neutrons, l_shield=l_wall, indices=indexes)
        self.signal_1 = curNuclei.calculate_signal(Nn=neutrons, Nph=photons, indices=indexes)
        self.signal_2 = curNuclei.calculate_signal(Nn=shielded_neutrons, Nph=shielded_photons, indices=indexes)

        self.answerLabel_1.setText(f'Полное число нейтронов без защиты: {curNuclei.answer_combing(self.signal_1)}')
        self.answerLabel_2.setText(f'Полное число нейтронов c защитой: {curNuclei.answer_combing(self.signal_2)}')

        self.spectrumEnergyGrid = curNuclei.FULL_ENERGY_GRID[indexes]
        self.photons = photons
        self.plot_spectrum(self.spectrumEnergyGrid, self.photons)
        self.photons = None
        #self.spectrumEnergyGrid = None

    # Plotters for signal and spectra
    def plot_signal(self):
        self.plotLayout.addWidget(self.signalCanvas)
        self.signalFigure.clear()
        self.signalCanvas.draw()
        self.ax = self.signalCanvas.figure.add_subplot(111)
        if self.signal_1 > 1:
            line1, = self.ax.plot(np.array([0, self.signal_1]), np.array([1.5, 1.5]), 'b',
                              label='Незащищенный детектор')
            self.ax.plot(np.array([0, 0]), np.array([0, 1.5]), 'b')
            self.ax.plot(np.array([self.signal_1, self.signal_1]), np.array([0, 1.5]), 'b')
        else:
            line1, = self.ax.plot(np.array([0, 0]), np.array([1.5, 1.5]), 'b',
                                  label='Незащищенный детектор (нет сигнала)')
            self.ax.plot(np.array([0, 0]), np.array([0, 1.5]), 'b')
            self.ax.plot(np.array([0, 0]), np.array([0, 1.5]), 'b')
        if self.signal_2 > 1:
            line2, = self.ax.plot(np.array([0, self.signal_2]), np.array([1.4, 1.4]), 'r',
                                      label='Защищенный детектор')
            self.ax.plot(np.array([0, 0]), np.array([0, 1.4]), 'r')
            self.ax.plot(np.array([self.signal_2, self.signal_2]), np.array([0, 1.4]), 'r')
        else:
            line2, = self.ax.plot(np.array([0, 0]), np.array([1.4, 1.4]), 'r',
                                  label='Защищенный детектор (нет сигнала)')
            self.ax.plot(np.array([0, 0]), np.array([0, 1.4]), 'r')
            self.ax.plot(np.array([0, 0]), np.array([0, 1.4]), 'r')
        self.ax.legend(handles=[line1, line2])
        self.ax.set_title(f'Сигнал на детекторах')
        self.signalCanvas.draw()
        # self.ax.set_xlabel(f'Энергия налетающего фотона, {self.EnergyModComboBox.currentText()}')
        # self.ax.set_ylabel('Сечение образования фотонейтрона, мбарн')        """

    def plot_spectrum(self, energy, photons, is_electron = False):
        if not self.spectrumPlotFlag:
            self.plotLayout.addWidget(self.spectrumCanvas)
            self.spectrumPlotFlag = True
        self.spectrumFigure.clear()
        self.spectrumCanvas.draw()
        self.ax_2 = self.spectrumCanvas.figure.add_subplot(111)
        self.ax_2.plot(energy, photons)
        self.ax_2.set_xlim(left=0, right=self.spectrumEnergyGrid[-1])
        self.ax_2.set_xlim(left=-5, right=20)
        self.ax_2.set_xlabel(f'Энергия фотонов, МэВ')
        self.ax_2.set_ylabel('Число фотонов, о.е.')
        self.ax_2.set_title(f'Спектр фотонов')
        self.spectrumCanvas.draw()


    # Manager of calculation, starting calculation depending on the method
    def Calculate(self):
        if self.methodComboBox.currentText() == "Расчёт от пучка электронов":
            self.beamCalculate()
        elif self.methodComboBox.currentText() == "Расчёт по спектру фотонов":
            self.spectrumCalculate()
        #print(f"Signal_1 = {self.signal_1}")
        #print(f"Signal_2 = {self.signal_2}")
        #self.plot_signal()

    # Test function
    def sensitivity_test(self):
        curNuclei = Nucleus(element="Fe", uplink=True)
        photons, indexes = curNuclei.preset_photons_spectrum(spectrum_type="TEST", k = 1, b = 1, c = 1, d = 1)
        neutrons = 10000
        self.signal_1 = curNuclei.calculate_signal(Nn=neutrons, Nph=photons, indices=indexes, is_electron=False)
        print(f'test signal = {self.signal_1}, energy = ')

    # A bunch of function controlling UI changes
    def filterChanged(self):
        if self.filterCheckBox.isChecked():
            self.aTextEdit.setVisible(True)
        else:
            self.aTextEdit.setVisible(False)

    def methodChange(self):
        for ph_widget, e_widget in zip(self.photon_widgets, self.electron_widgets):
            if isinstance(ph_widget, (QTextEdit)):
                ph_widget.setText('')
            if isinstance(e_widget, (QTextEdit)):
                e_widget.setText('')

        if self.methodComboBox.currentText() == "Расчёт от пучка электронов":
            for widget in self.photon_widgets:
                widget.setVisible(False)
            for widget in self.electron_widgets:
                widget.setVisible(True)

        elif self.methodComboBox.currentText() == "Расчёт по спектру фотонов":
            for widget in self.photon_widgets:
                widget.setVisible(True)
            for widget in self.electron_widgets:
                widget.setVisible(False)
            self.spectrumComboBox.setCurrentText("Линейный")
            self.spectrumChange()

    def recalculateGridCheckstateChanged(self):
        if self.recalculateCheckBox.isChecked():
            self.minTextEdit.setVisible(True)
            self.maxTextEdit.setVisible(True)
            self.stepTextEdit.setVisible(True)

            self.minLabel.setVisible(True)
            self.maxLabel.setVisible(True)
            self.stepLabel.setVisible(True)

            self.recalculateFlag = True
        else:
            self.minTextEdit.setVisible(False)
            self.maxTextEdit.setVisible(False)
            self.stepTextEdit.setVisible(False)

            self.minLabel.setVisible(False)
            self.maxLabel.setVisible(False)
            self.stepLabel.setVisible(False)

            self.recalculateFlag = False

    def spectrumChange(self):
        if self.spectrumComboBox.currentText() == "Линейный":
            self.cLabel.setVisible(False)
            self.dLabel.setVisible(False)
            self.cTextEdit.setVisible(False)
            self.dTextEdit.setVisible(False)
            self.equationLabel.setText("N = aE + b")
        elif self.spectrumComboBox.currentText() == "Параболический":
            self.cLabel.setVisible(True)
            self.dLabel.setVisible(False)
            self.cTextEdit.setVisible(True)
            self.dTextEdit.setVisible(False)
            self.cTextEdit.clear()
            self.equationLabel.setText("N = aE^2 + bE + c")
        elif self.spectrumComboBox.currentText() == "Гаусс":
            self.cLabel.setVisible(True)
            self.dLabel.setVisible(True)
            self.cTextEdit.setVisible(True)
            self.dTextEdit.setVisible(True)
            self.cTextEdit.clear()
            self.equationLabel.setText("N = ae^((E-b)^2/c^2) + d")

    def valueFromText(self, string):
        l_string = string.replace(",", ".")
        l_string = l_string.replace(" ", "")
        str_list = l_string.split("e")
        if not l_string:
            return 0
        if l_string[0].isnumeric():
            try:
                number = float(float(str_list[0]) * (10 ** float(str_list[1])))
            except IndexError:
                number = float(str_list[0])
            return number
        elif l_string[0] == "e":
            number = 10 ** float(str_list[1])
            return number
        elif l_string[0] == "-":
            if l_string[1].isnumeric():
                try:
                    number = float(float(str_list[0]) * (10 ** float(str_list[1])))
                except IndexError:
                    number = float(str_list[0])
                return number
            elif l_string[1] == "e":
                number = -10 ** float(str_list[1])
                return number
        else:
            self.errorPopup = 1

    # Function controlling plot actions
    def plot_connect(self):
        self.cidpress = self.spectrumFigure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.spectrumFigure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.spectrumFigure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.ax_2.axes:
            return
        self.press = event.xdata, event.ydata

    def on_motion(self, event):
        if self.press is None or event.inaxes != self.ax_2.axes:
            return
        if not self.patchFlag:
            self.patchFlag = True
            if event.xdata > self.press[0]:
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax_2.get_ylim()[1], color=('blue', 0.3))
                self.ax_2.add_patch(self.rectangle)
                self.spectrumCanvas.draw()
            else:
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax_2.get_ylim()[1], color=('red', 0.4))
                self.ax_2.add_patch(self.rectangle)
                self.spectrumCanvas.draw()
        else:
            if event.xdata > self.press[0]:
                self.rectangle.remove()
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax_2.get_ylim()[1],
                                           color=('blue', 0.3))
                self.ax_2.add_patch(self.rectangle)
                self.spectrumCanvas.draw()
            else:
                self.rectangle.remove()
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax_2.get_ylim()[1],
                                           color=('red', 0.4))
                self.ax_2.add_patch(self.rectangle)
                self.spectrumCanvas.draw()

    def on_release(self, event):
        if event.inaxes != self.ax_2.axes:
            return
        if self.press[0] < event.xdata:
            self.rectangle.remove()
            #self.plot_spectrum(self.spectrumEnergyGrid, self.photons)
            self.ax_2.set_xlim(left=self.press[0], right=event.xdata)
            self.spectrumCanvas.draw()
        else:
            self.rectangle.remove()
            self.ax_2.set_xlim(left=0, right= 10**8 + 100000)
            self.spectrumCanvas.draw()
            #self.spectrumFigure.clear()
            #self.plot_spectrum(self.spectrumEnergyGrid, self.photons)
            #self.plot_signal()
        self.press = None
        self.patchFlag = False

    def disconnect(self):
        self.spectrumFigure.canvas.mpl_disconnect(self.cidpress)
        self.spectrumFigure.canvas.mpl_disconnect(self.cidrelease)
        self.spectrumFigure.canvas.mpl_disconnect(self.cidmotion)
