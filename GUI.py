import sys, os
from dbm import error

import numpy as np
import xlsxwriter as xl
from PyQt6.QtWidgets import QWidget, QMainWindow, QDialog, QDialogButtonBox, QVBoxLayout, QLabel, QApplication, \
    QFileDialog, QStackedWidget, QCheckBox, QComboBox, QPushButton
from PyQt6 import uic
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from numpy.f2py.auxfuncs import isstring
from Counter import *
from nuclei import *


energies = {
    "МэВ" : pow(10,6),
    "кэВ": pow(10,3),
    "эВ": 1
}

progs = {
    'Data reading' : 0,
    'Photoneutron calculations' : 1,
    'Journal notes' : 2
}

spectra = {
    "Линейный",
    "Параболический",
    "Гаусс"
}

ui_photon, _ = uic.loadUiType('Elphoneut.ui')
ui_counter, _ = uic.loadUiType('Neutrongit')
ui_window, _ = uic.loadUiType('Journal.ui')

class JournalWindow(QMainWindow, ui_window):
    def __init__(self, parent=None):
        super(JournalWindow, self).__init__(parent)
        self.setupUi(self)
        self.press = None
        self.save = None
        self.filename = None
        self.savename = 'NeutronData.xlsx'
        self.checkboxFlag = 0
        self.invertcoeff = 1


        self.invertcheckBox.setVisible(False)
        self.plotFigure = Figure()
        self.Canvas = FigureCanvas(self.plotFigure)
        self.verticalLayout_2.addWidget(self.Canvas)

        self.ThresholdSpinBox_1.setValue(0.01)
        self.ThresholdSpinBox_2.setValue(0.01)

        for i in progs.keys():
            self.ProgramChannelBox.addItem(i)
        self.ProgramChannelBox.setCurrentText('Journal notes')

        self.CalculateButton.clicked.connect(self.calculate)
        self.impulseEdit.textEdited.connect(self.newfilereset)


        self.label_I.setText(f'')
        self.neutronEdit.setPlaceholderText('Путь к файлу XR')
        self.currentEdit.setPlaceholderText('Путь к файлу с токами')
        self.gyroEdit.setPlaceholderText('Путь к файлу с гиротронами')

    def newfilereset(self):
        self.PlotFlag = 0
        self.plotFigure.clear()
        self.Canvas.draw()


    def calculate(self):
        if ERROR == 1:
            self.errorPopup('Cyrillic')
            return

        self.label_SXR1.setText(f'')
        self.label_SXR2.setText(f'')
        self.label_SXR3.setText(f'')
        self.label_SXR4.setText(f'')
        self.label_HXR1.setText(f'')
        self.label_I.setText(f'')

        if self.impulseEdit.text() == '':
            self.errorPopup('NoImpulse')
            return

        self.number = self.impulseEdit.text()
        if self.neutronEdit.text() == '':
            self.neutfile = params['journalData']['paths']['1'] + '/'  + params['journalData']['filenames']['1'].partition('number')[0]  + self.number + params['journalData']['filenames']['1'].partition('number')[2]
        else:
            self.neutfile = (self.neutronEdit.text() + '/'  + params['journalData']['filenames']['1'].partition('number')[0] + self.number + params['journalData']['filenames']['1'].partition('number')[2])
        if self.currentEdit.text() == '':
            self.Ifile = params['journalData']['paths']['2'] + '/'  +  self.number + params['journalData']['filenames']['2']
        else:
            self.Ifile = (self.currentEdit.text() + '/'   +  self.number + params['journalData']['filenames']['2'])
        try:
            self.neutCalc([self.neutfile])
            self.peakCalc(self.neutfile,'SXR_Mon_')
            self.peakCalc(self.neutfile, 'HXR_Mon_')
        except FileNotFoundError:
            self.errorPopup('XR')
        try:
            self.peakCalc(self.Ifile, 'ipl')
        except FileNotFoundError:
            self.errorPopup('Ipl')
        try:
            self.tdms_file = td.TdmsFile.read(self.neutfile)
            self.toggles = np.zeros(len(self.tdms_file['DAQmx'].channels()) - 1)
            self.invertcheckBox.setVisible(True)
            self.invertcheckBox.checkStateChanged.connect(self.plotinvert)
            if not self.checkboxFlag:
                self.checkboxes = []
                for group in self.tdms_file['DAQmx'].channels():
                    if group.name != '0.inf':
                        checkbox = QCheckBox(f"{group.name}")
                        self.checkboxes.append(checkbox)
                        self.verticalLayout.addWidget(checkbox)
                        self.checkboxFlag = 1
                    checkbox.checkStateChanged.connect(self.checkTick)
        except FileNotFoundError:
            self.invertcheckBox.setVisible(False)
            self.errorPopup('XR')


    def plotinvert(self):
        try:
            self.invertcoeff = self.invertcoeff * -1
            self.tickDraw()
        except FileNotFoundError:
            self.errorPopup('XR')

    def checkTick(self):
        i = 0
        for check in self.checkboxes:
            if check.isChecked():
                self.toggles[i] = 1
            else:
                self.toggles[i] = 0
            i += 1
        self.tickDraw()

    def tickDraw(self):
        if self.PlotFlag:
            self.xLimit = self.ax.get_xlim()
        self.plotFigure.clear()
        self.Canvas.draw()
        self.ax = self.Canvas.figure.add_subplot(111)
        try:
            for i in range(0,len(self.toggles)):
                if self.toggles[i]:
                    count = Calculator(self.neutfile,noise_zero=self.ThresholdSpinBox_1.value(),noise_peak=self.ThresholdSpinBox_2.value(), channel=self.tdms_file['DAQmx'].channels()[i].name)
                    self.data = count.true_data*self.invertcoeff
                    self.ax.plot(self.data)
                    self.data = None
                    if self.tdms_file['DAQmx'].channels()[i].name == 'Neut_Mon_1' or self.tdms_file['DAQmx'].channels()[i].name == 'Neut_Mon_2':
                        self.peaks, self.data_peaks, self.amount = count.find_peaks()
                        self.ax.plot(self.peaks, self.data_peaks*self.invertcoeff, "x")
                        self.peaks,self.data_peaks, self.amount  = None, None, None
        except FileNotFoundError:
            self.errorPopup('XR')
        if self.PlotFlag:
            self.ax.set_xlim(self.xLimit)
        if not self.PlotFlag:
            self.PlotFlag = 1
        self.Canvas.draw()
        self.connect()

    def peakCalc(self, file, channel):
        if channel == 'SXR_Mon_':
            signalSXR = np.zeros(4)
            for i in range(0,len(signalSXR)):
                signal = Signaler(file, 'MIV', channel + str(i+1))
                signalSXR[i] = signal.peak()
            self.label_SXR1.setText(f'{signalSXR[0]} мкВт/см^2')
            self.label_SXR2.setText(f'{signalSXR[1]} мкВт/см^2')
            self.label_SXR3.setText(f'{signalSXR[2]} мкВт/см^2')
            self.label_SXR4.setText(f'{signalSXR[3]} мкВт/см^2')

        elif channel == 'HXR_Mon_':
            signalHXR = np.zeros(1)
            for i in range(0, len(signalHXR)):
                signal = Signaler(file, 'MIV', channel + str(i + 1))
                signalHXR[i] = signal.peak()
                self.label_HXR1.setText(f'{signalHXR[0]} мкВт')
        elif channel == 'ipl':
            signal = Signaler(file, 'EMD', channel)
            signalIpl = signal.peak()
            self.label_I.setText(f'{signalIpl} кА')

    def neutCalc(self, file):
        inCalc = MultiCalcultor(file, noise_zero=self.ThresholdSpinBox_1.value(),
                                 noise_peak=self.ThresholdSpinBox_2.value())
        amount = inCalc.find_peaks()
        self.label_Neut1.setText(f'Mon_1: {amount[0]}')
        self.label_Neut2.setText(f'Mon_2: {amount[1]}')

    def errorPopup(self, error):
        dlg = QDialog(self)
        dlg.setWindowTitle("Help")
        buttonBox = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        buttonBox.clicked.connect(dlg.close)
        layout = QVBoxLayout()
        if error == 'XR':
            layout.addWidget(
                QLabel(f"В папке нет файла XR с номером {self.number}."))
        if error == 'Ipl':
            layout.addWidget(
                QLabel(f"В папке нет файла с токами с номером {self.number}."))
        if error == 'NoImpulse':
            layout.addWidget(
                QLabel(f"Импульс не выбран"))
        if error == 'Cyrillic':
            layout.addWidget(
                QLabel(f"На пути к файлу есть Кириллица"))
        layout.addWidget(buttonBox)
        dlg.setLayout(layout)
        dlg.exec()


    def connect(self):
        self.cidpress = self.plotFigure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.plotFigure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.plotFigure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.ax.axes:
            return
        self.press = event.xdata, event.ydata
        self.temp_line = []
        self.lim = self.ax.get_ylim()
        self.ax.plot([self.press[0], self.press[0]], [self.lim[0], self.lim[1]])
        self.Canvas.draw()

    def on_motion(self, event):
        if self.press is None or event.inaxes != self.ax.axes:
            return
        #       to comment
        #self.ax.clear()
        #self.plot_graph()
        #self.ax.plot([self.press[0], self.press[0]], [0, self.lim])
        #self.temp = event.xdata, event.ydata
        #self.temp_line = self.ax.plot([self.temp[0], self.temp[0]], [0, self.lim])
        #self.Canvas.draw()
        #       stop comment

    def on_release(self, event):
        if event.inaxes != self.ax.axes:
            return
        if self.press[0] < event.xdata:
            #       to comment
            self.plotFigure.clear()
            self.tickDraw()
            #       stop comment
            self.ax.set_xlim(left=self.press[0], right=event.xdata)
        else:
            self.PlotFlag = 0
            #self.ax.set_xlim(left=0, right=len(self.data))
            self.plotFigure.clear()
            self.tickDraw()
        self.press = None
        self.Canvas.draw()

    def disconnect(self):
        self.plotFigure.canvas.mpl_disconnect(self.cidpress)
        self.plotFigure.canvas.mpl_disconnect(self.cidrelease)
        self.plotFigure.canvas.mpl_disconnect(self.cidmotion)

class CounterWindow(QMainWindow, ui_counter):
    def __init__(self, parent=None):
        super(CounterWindow, self).__init__(parent)
        self.setupUi(self)
        self.patchFlag = False
        self.press = None
        self.save = None
        self.filename = None
        self.savename = 'NeutronData.xlsx'
        self.invertCoef = 1

        self.ThresholdSpinBox_1.setValue(0.01)
        self.ThresholdSpinBox_2.setValue(0.01)

        for i in progs.keys():
            self.ProgramChannelBox.addItem(i)
        self.ProgramChannelBox.setCurrentText('Data reading')

        self.plotFigure = Figure()
        self.plot_connect()
        self.Canvas = FigureCanvas(self.plotFigure)
        self.verticalLayout.addWidget(self.Canvas)

        self.OpenButton.clicked.connect(self.openFile)
        self.HelpButton.clicked.connect(self.helpPopup)
        self.addButton.clicked.connect(self.addChannel)
        self.pleaseButton.clicked.connect(self.please)


        self.lineEdit.setPlaceholderText('Название файла')

        self.calculateButton.setVisible(False)
        self.invertCheckBox.setVisible(False)
        self.integrateCheckBox.setVisible(False)
        self.hideCheckBox.setVisible(False)
        self.peaksCheckBox.setVisible(False)
        self.normCheckBox.setVisible(False)
        self.addButton.setVisible(False)
        self.ChannelComboBox.setVisible(False)

    def addChannel(self):
        newButton = QPushButton('Убрать канал')
        newComboBox = QComboBox()
        newCheckBox = QCheckBox('Инвертировать')
        newCheckBox_2 = QCheckBox("Интегрировать")
        hideCheckBox = QCheckBox('Скрыть')
        normCheckBox = QCheckBox('Нормировать')
        newButton.clicked.connect(self.create_removeChannel(newButton, newComboBox, newCheckBox, newCheckBox_2,
                                                            hideCheckBox, normCheckBox))
        for group in self.tdms_file.groups():
            for channel in group.channels():
                newComboBox.addItem(channel.name)
                print(channel.name)

        self.removeButtonLayout.addWidget(newButton, alignment=Qt.AlignmentFlag.AlignTop)
        self.comboBoxLayout.addWidget(newComboBox, alignment=Qt.AlignmentFlag.AlignTop)
        self.checkBoxLayout.addWidget(newCheckBox, alignment=Qt.AlignmentFlag.AlignTop)
        self.checkBoxLayout_2.addWidget(newCheckBox_2, alignment=Qt.AlignmentFlag.AlignTop)
        self.hideCheckBoxLayout.addWidget(hideCheckBox, alignment=Qt.AlignmentFlag.AlignTop)
        self.normCheckBoxLayout.addWidget(normCheckBox, alignment=Qt.AlignmentFlag.AlignTop)

    def create_removeChannel(self, button, combobox, checkbox, checkbox_2, checkbox_3, checkbox_4):
        def removeChannel():
            self.removeButtonLayout.removeWidget(button)
            self.comboBoxLayout.removeWidget(combobox)
            self.checkBoxLayout.removeWidget(checkbox)
            self.checkBoxLayout_2.removeWidget(checkbox_2)
            self.hideCheckBoxLayout.removeWidget(checkbox_3)
            self.normCheckBoxLayout.removeWidget(checkbox_4)
            button.setVisible(False)
            combobox.setVisible(False)
            checkbox.setVisible(False)
            checkbox_2.setVisible(False)
            checkbox_3.setVisible(False)
            checkbox_4.setVisible(False)
        return removeChannel

    def please(self):
        inCalc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                            noise_peak=self.ThresholdSpinBox_2.value(),
                            channel="Neut_Mon_1",
                            inverted=-1)
        peaks, data_peaks, amount_1 = inCalc.find_peaks()
        inCalc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                            noise_peak=self.ThresholdSpinBox_2.value(),
                            channel="Neut_Mon_2",
                            inverted=-1,
                            diff=True)
        peaks, data_peaks, amount_2 = inCalc.find_peaks()
        if amount_2 != 0:
            self.label_2.setText(f'{amount_1 / amount_2}')

    def openFile(self):
        self.dialog = QFileDialog(self)
        self.filename, _ = self.dialog.getOpenFileNames(self, 'Open file', './', '(*.tdms)')
        if self.filename == '' or self.filename is None or self.filename == []:
            return
        if len(self.filename) != 1:
            #for item in self.ChannelComboBox.items():
            #    self.ChannelComboBox.removeItem(item)
            #self.ChannelComboBox.setCurrentText('')
            self.OpenButton.setText('Выбрано несколько файлов')
            self.calculateButton.setText('Рассчитать нейтроны в каждом файле')
            self.calculateButton.clicked.connect(self.dataWrite)
            self.calculateButton.setVisible(True)
            #self.label_2.setText(f'')
            self.multifile = True
            #self.ChannelComboBox.setVisible(False)
        else:
            #for i in channels:
            #    self.ChannelComboBox.addItem(i)
            self.filename = self.filename[0]
            self.OpenButton.setText('Текущий файл:' + self.filename.split('/')[len(self.filename.split('/'))-1])
            self.calculateButton.setText('Построить график')
            self.tdms_file = td.TdmsFile.read(self.filename)
            for group in self.tdms_file['Neut'].channels():
                self.ChannelComboBox.addItem(group.name)
            self.multifile = False
            if self.filename == None or self.filename == []:
                print(f'Ошибка: Файл не выбран')
                return
            self.ChannelComboBox.setVisible(True)
            self.invertCheckBox.setVisible(True)
            self.calculateButton.setVisible(True)
            self.integrateCheckBox.setVisible(True)
            self.hideCheckBox.setVisible(True)
            self.peaksCheckBox.setVisible(True)
            self.normCheckBox.setVisible(True)
            self.addButton.setVisible(True)
            self.calculateButton.clicked.connect(self.plot_simple)
            inCalc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                                     noise_peak=self.ThresholdSpinBox_2.value(),
                                     channel="Neut_Mon_1",
                                     inverted=-1)
            peaks, data_peaks, amount_1 = inCalc.find_peaks()
            inCalc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                                     noise_peak=self.ThresholdSpinBox_2.value(),
                                     channel="Neut_Mon_2",
                                     inverted=-1)
            peaks, data_peaks, amount_2 = inCalc.find_peaks()

            val_1 = 666666
            if amount_2 != 0:
                val_1 = amount_1 / amount_2

            inCalc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                                     noise_peak=self.ThresholdSpinBox_2.value(),
                                     channel="Neut_Mon_1",
                                     inverted=1)
            peaks, data_peaks, amount_1 = inCalc.find_peaks()
            inCalc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                                     noise_peak=self.ThresholdSpinBox_2.value(),
                                     channel="Neut_Mon_2",
                                     inverted=1)
            peaks, data_peaks, amount_2 = inCalc.find_peaks()
            val_2 = 666666
            if amount_2 != 0:
                val_2 = amount_1 / amount_2

            self.label_2.setText(f'inv:{round(val_1, 2)}, non_inv:{round(val_2, 2)}')
        #self.plotFigure.clear()
        #self.Canvas.draw()

    def invertCheckBoxCoefficient(self, checkbox):
        if checkbox.isChecked():
            return -1
        else:
            return 1

    def stateChanged(self):
        self.inCalc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                                 noise_peak=self.ThresholdSpinBox_2.value(),
                                 channel=self.ChannelComboBox.currentText(),
                                 inverted=self.invertCheckBoxCoefficient(self.invertCheckBox))
    def helpPopup(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("Help")
        buttonBox = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        buttonBox.clicked.connect(dlg.close)
        layout = QVBoxLayout()
        layout.addWidget(QLabel(
            "В правом меню можно переключаться между считыванием данных нейтронных датчиков и расчётом фотонейтронов."))
        layout.addWidget(QLabel("В поле Порог от нуля устанавливается минимальная высота пика от нуля."))
        layout.addWidget(QLabel("В поле Порог на пике устанавливается минимальная высота пика от соседних значений."))
        layout.addWidget(
            QLabel("Перед запуском программы необходимо открыть .tdlm файл или нажав на соответсвующую кнопку."))
        layout.addWidget(
            QLabel("При выборе нескольких файлов будет рассчитано число нейтронов в каждом канале."))
        layout.addWidget(
            QLabel("В таком случае число нейтронов в каждом канале будет сохранено в файл с названием из соответствующего поля."))
        layout.addWidget(
            QLabel(
                "Если название не будет указано - файл сохранится под стандартным названием NeutronData.xlsx"))
        layout.addWidget(
            QLabel("При выборе одного файла будет построен график имеющегося сигнала."))
        layout.addWidget(QLabel("Зажав мышку на графике и проведя слева направо можно увеличить выделенный интервал."))
        layout.addWidget(
            QLabel("Зажав мышку на графике и проведя справа налево можно вернуть интервал к стандартному."))

        layout.addWidget(buttonBox)
        dlg.setLayout(layout)
        dlg.exec()
    def plotter(self):
        return

    def decimate_minmax(self, x, y, bins=1000):
        if len(x) <= bins:
            return x, y

        indices = np.linspace(0, len(x) - 1, bins + 1, dtype=int)
        new_x, new_y = [], []

        for i in range(len(indices) - 1):
            segment = slice(indices[i], indices[i + 1])
            x_seg, y_seg = x[segment], y[segment]
            min_idx = np.argmin(y_seg) + indices[i]
            max_idx = np.argmax(y_seg) + indices[i]
            new_x.extend([x[min_idx], x[max_idx]])
            new_y.extend([y[min_idx], y[max_idx]])

        return np.array(new_x), np.array(new_y)


    def plot_simple(self):
        self.stateChanged()
        self.plotFigure.clear()
        self.Canvas.draw()
        handles = []
        self.ax = self.Canvas.figure.add_subplot(111)
        if not self.hideCheckBox.isChecked():
            if self.integrateCheckBox.isChecked():
                vals = self.inCalc.integrate(self.ThresholdSpinBox_3.value())
                if self.normCheckBox.isChecked():
                    maxim = max(vals)
                    for i in range(0, len(vals)):
                        vals[i] = vals[i]/abs(maxim)
            elif self.differentiateCheckBox.isChecked():
                vals = self.inCalc.differentiate()
                if self.normCheckBox.isChecked():
                    maxim = max(vals)
                    for i in range(0, len(vals)):
                        vals[i] = vals[i] / abs(maxim)
            else:
                print(5)
                vals = self.inCalc.true_data()
                if self.normCheckBox.isChecked():
                    maxim = max(vals)
                    for i in range(0, len(vals)):
                        vals[i] = vals[i]/abs(maxim)

            line, = self.ax.plot(self.inCalc.time, vals, 'b', label=self.ChannelComboBox.currentText())
            if self.peaksCheckBox.isChecked():
                peaks, data_peaks, amount = self.inCalc.find_peaks()
                self.ax.plot(self.inCalc.time[peaks], data_peaks, "x")
            handles.append(line)

        for i in range(0, self.removeButtonLayout.count()):
            if not self.hideCheckBoxLayout.itemAt(i).widget().isChecked():
                channelName = self.comboBoxLayout.itemAt(i).widget().currentText()
                invertCheck = self.checkBoxLayout.itemAt(i).widget()
                Calc = Calculator(self.filename, noise_zero=self.ThresholdSpinBox_1.value(),
                                     noise_peak=self.ThresholdSpinBox_2.value(),
                                     channel=channelName,
                                     inverted=self.invertCheckBoxCoefficient(invertCheck))
                if self.checkBoxLayout_2.itemAt(i).widget().isChecked():
                    chanData = Calc.integrate(self.ThresholdSpinBox_3.value())
                else:
                    chanData = Calc.true_data()
                if self.normCheckBoxLayout.itemAt(i).widget().isChecked():
                    maxim = max(chanData)
                    for j in range(0, len(chanData)):
                        chanData[j] = chanData[j] / abs(maxim)
                if self.comboBoxLayout.itemAt(i).widget().currentText() == "Neut_Mon_2":
                    line, = self.ax.plot(self.inCalc.time, chanData, 'r',label=channelName)
                elif self.comboBoxLayout.itemAt(i).widget().currentText() == "SXR_Mon_4":
                    line, = self.ax.plot(self.inCalc.time, chanData, 'g', label=channelName)
                elif self.comboBoxLayout.itemAt(i).widget().currentText() == "HXR_Mon_1":
                    line, = self.ax.plot(self.inCalc.time, chanData, 'm', label=channelName)
                else:
                    line, = self.ax.plot(self.inCalc.time, chanData, label=channelName)
                handles.append(line)

        self.ax.legend(handles=handles)
        self.Canvas.draw()
        self.plot_connect()

    def dataWrite(self):
        print(1)
        if self.savename == '':
            self.savename = 'CalculatedData.xlsx'
        print(2)
        output = xl.Workbook(self.savename)
        outputSheet = output.add_worksheet()
        index = 0
        print(3)
        for i in range(len(self.filename)):
            file = td.TdmsFile.read(self.filname[i])
            monitors = True
            outputSheet.write(0, index*3+1, self.filename[i])
            all_channels = []
            for group in file.groups():
                for chan in group.channels():
                    all_channels.append(chan.name)
            while monitors:
                neut_name = 'Neut_mon_'
                i = 1
                if (neut_name + str(i)) in all_channels:
                    outputSheet.write(1, index+i, neut_name + str(i))
                    inCalc = Calculator(file, noise_zero=self.ThresholdSpinBox_1.value(),
                                    noise_peak=self.ThresholdSpinBox_2.value(),
                                    channel=neut_name + str(i),
                                    inverted=1)
                    peaks, data_peaks, amount = inCalc.find_peaks()
                    outputSheet.write(2, index + i, amount)
                else:
                    monitors = False
                i += 1

        outputSheet.autofit()
        outputSheet.center_vertically()
        output.close()

    def plot_connect(self):
        self.cidpress = self.plotFigure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.plotFigure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.plotFigure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.ax.axes:
            return
        self.press = event.xdata, event.ydata

    def on_motion(self, event):
        if self.press is None or event.inaxes != self.ax.axes:
            return
        if not self.patchFlag:
            self.patchFlag = True
            if event.xdata > self.press[0]:
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax.get_ylim()[1],
                                           color=('blue', 0.3))
                self.ax.add_patch(self.rectangle)
                self.Canvas.draw()
            else:
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax.get_ylim()[1],
                                           color=('red', 0.4))
                self.ax.add_patch(self.rectangle)
                self.Canvas.draw()
        else:
            if event.xdata > self.press[0]:
                self.rectangle.remove()
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax.get_ylim()[1],
                                           color=('blue', 0.3))
                self.ax.add_patch(self.rectangle)
                self.Canvas.draw()
            else:
                self.rectangle.remove()
                self.rectangle = Rectangle((self.press[0], 0), event.xdata - self.press[0], self.ax.get_ylim()[1],
                                           color=('red', 0.4))
                self.ax.add_patch(self.rectangle)
                self.Canvas.draw()

    def on_release(self, event):
        if event.inaxes != self.ax.axes:
            return
        if self.press[0] < event.xdata:
            self.rectangle.remove()
            self.ax.set_xlim(left=self.press[0], right=event.xdata)
            self.Canvas.draw()
        else:
            self.rectangle.remove()
            self.plotFigure.clear()
            self.Canvas.draw()
            self.plot_simple()
            self.stateChanged()
            # self.spectrumFigure.clear()
            # self.plot_spectrum(self.spectrumEnergyGrid, self.photons)
            # self.plot_signal()
        self.press = None
        self.patchFlag = False

    def plot_disconnect(self):
        self.plotFigure.canvas.mpl_disconnect(self.cidpress)
        self.plotFigure.canvas.mpl_disconnect(self.cidrelease)
        self.plotFigure.canvas.mpl_disconnect(self.cidmotion)

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

        self.methodComboBox.addItem('Расчёт от пучка электронов')
        self.methodComboBox.addItem('Расчёт по спектру фотонов')
        self.methodComboBox.setCurrentText('Расчёт от пучка электронов')
        self.methodChange()

        #self.helpButton.clicked.connect(self.helpPopup)

        for spectrum in spectra:
            self.spectrumComboBox.addItem(spectrum)
        for target in data.keys():
            self.targetComboBox.addItem(target)
        self.targetComboBox.setCurrentText("Fe")

        self.recalculateGridCheckstateChanged()

        self.signalFigure = Figure()
        self.signalCanvas = FigureCanvas(self.signalFigure)
        self.signalCanvas.setMaximumHeight(400)
        self.signalCanvas.setMinimumHeight(300)
        self.signalCanvas.setMinimumWidth(450)

        self.spectrumFigure = Figure()
        self.spectrumCanvas = FigureCanvas(self.spectrumFigure)
        self.spectrumCanvas.setMaximumHeight(400)
        self.spectrumCanvas.setMinimumHeight(300)
        self.spectrumCanvas.setMinimumWidth(450)

        self.plot_connect()

        self.recalculateCheckBox.stateChanged.connect(self.recalculateGridCheckstateChanged)
        self.filterCheckBox.stateChanged.connect(self.filterChanged)

        self.methodComboBox.currentTextChanged.connect(self.methodChange)
        self.spectrumComboBox.currentTextChanged.connect(self.spectrumChange)


        self.calculatePushButton.clicked.connect(self.Calculate)

    def filterChanged(self):
        if self.filterCheckBox.isChecked():
            self.aTextEdit.setVisible(True)
        else:
            self.aTextEdit.setVisible(False)

    def sensitivity_test(self):
        curNuclei = Nucleus(element="Fe", uplink=True)
        photons, indexes = curNuclei.preset_photons_spectrum(spectrum_type="TEST", k = 1, b = 1, c = 1, d = 1)
        neutrons = 10000
        self.signal_1 = curNuclei.calculate_signal(Nn=neutrons, Nph=photons, indices=indexes, is_electron=False)
        print(f'test signal = {self.signal_1}, energy = ')

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

    def methodChange(self):
        if self.methodComboBox.currentText() == "Расчёт от пучка электронов":
            self.enableElectrons()
        elif self.methodComboBox.currentText() == "Расчёт по спектру фотонов":
            self.enablePhotons()

    def enableElectrons(self):
        self.aLabel.setVisible(False)
        self.bLabel.setVisible(False)
        self.cLabel.setVisible(False)
        self.dLabel.setVisible(False)
        self.aTextEdit.setVisible(False)
        self.bTextEdit.setVisible(False)
        self.cTextEdit.setVisible(False)
        self.dTextEdit.setVisible(False)
        self.spectrumComboBox.setVisible(False)
        self.equationLabel.setText('')
        self.equationLabel.setVisible(False)

        self.filterCheckBox.setVisible(True)
        self.energyLabel.setVisible(True)
        self.amountLabel.setVisible(True)
        self.amountTextEdit.setVisible(True)
        self.energyTextEdit.setVisible(True)

    def enablePhotons(self):
        self.spectrumComboBox.setVisible(True)
        self.equationLabel.setVisible(True)
        self.aTextEdit.clear()
        self.bTextEdit.clear()
        self.cTextEdit.clear()
        self.dTextEdit.clear()
        self.spectrumComboBox.setCurrentText("Линейный")
        self.spectrumChange()
        self.aLabel.setVisible(True)
        self.bLabel.setVisible(True)
        self.aTextEdit.setVisible(True)
        self.bTextEdit.setVisible(True)
        self.filterCheckBox.setVisible(False)
        self.energyLabel.setVisible(False)
        self.amountLabel.setVisible(False)
        self.amountTextEdit.setVisible(False)
        self.energyTextEdit.setVisible(False)

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
            curNuclei = Nucleus(E0=E0, Ne=Ne, element=self.curElement, uplink=False)
            curNuclei.recalculate_crosssections()
            curNuclei = Nucleus(E0=E0, Ne=Ne, element=curElement, uplink=True, is_electron=True)
        else:
            curNuclei = Nucleus(E0=E0, Ne=Ne, element=curElement, uplink=True, is_electron=True)

        if self.targetComboBox.currentText() != 'Сталь 12Х18Е10Т':
            print(self.valueFromText(self.aTextEdit.toPlainText()))
            self.photons = curNuclei.calculate_photons(filtering=self.valueFromText(self.aTextEdit.toPlainText()))
            self.spectrumEnergyGrid = curNuclei.PhotonEnergies
            self.plot_spectrum(self.spectrumEnergyGrid, self.photons, is_electron=True)
            print(7)
            neutrons = curNuclei.calculate_neutrons(Nph=self.photons, l_wall=l)
            #R = 1.5
            #S = 0.5*0.15
            #neutrons = neutrons*S/(4*3.14159*R*R)
            print(8)
            shielded_photons = curNuclei.calculate_shielded_photons(Nph=self.photons, l_shield=l_wall)
            shielded_neutrons = curNuclei.calculate_shielded_neutrons(Nn=neutrons, l_shield=l_wall)
            print(9)
            #shielded_neutrons = shielded_neutrons * S / (4 * 3.14159 * R * R)
            self.signal_1 = curNuclei.calculate_signal(Nn=neutrons, Nph=self.photons)
            self.signal_2 = curNuclei.calculate_signal(Nn=shielded_neutrons, Nph=shielded_photons)
            print(10)
            #self.answerLabel_1.setText(f'Полное число сигналов без защиты: {curNuclei.answer_combing(self.signal_1)}')
            #self.answerLabel_2.setText(f'Полное число сигналов c защитой: {curNuclei.answer_combing(self.signal_2)}')
            #self.answerLabel_3.setText(f'Соотношение числа сигналов без/c защитой: {round(self.signal_1/self.signal_2, 3)}')
        else:
            neutrons = curNuclei.answer_combing(curNuclei.alloy_calculate_photoneutrons(Ne=Ne, E0=E0, l=l))
            self.ResultLabel.setText(f'Полное число фотонейтронов: {neutrons}')

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

    def Calculate(self):
        if self.methodComboBox.currentText() == "Расчёт от пучка электронов":
            self.beamCalculate()
        elif self.methodComboBox.currentText() == "Расчёт по спектру фотонов":
            self.spectrumCalculate()
        #print(f"Signal_1 = {self.signal_1}")
        #print(f"Signal_2 = {self.signal_2}")
        #self.plot_signal()

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
        loc_en = energy*10**-6
        loc_en[-8] = 16.25
        loc_en[-7] = 16.5
        loc_en[-6] = 16.75
        loc_en[-5] = 17
        loc_en[-4] = 17.25
        loc_en[-3] = 17.5
        loc_en[-2] = 17.75
        loc_en[-1] = 18
        loc_phot = photons
        last = loc_phot[-1]
        loc_phot[-8] = last*0.975
        loc_phot[-7] = last*0.925
        loc_phot[-6] = last*0.850
        loc_phot[-5] = last*0.750
        loc_phot[-4] = last*0.600
        loc_phot[-3] = last*0.350
        loc_phot[-2] = last*0.100
        loc_phot[-1] = 0
        print(1)
        if not self.spectrumPlotFlag:
            self.plotLayout.addWidget(self.spectrumCanvas)
            self.spectrumPlotFlag = True
        self.spectrumFigure.clear()
        self.spectrumCanvas.draw()
        print(2)
        self.ax_2 = self.spectrumCanvas.figure.add_subplot(111)
        self.ax_2.plot(loc_en, loc_phot)
        self.ax_2.set_xlim(left=0, right=self.spectrumEnergyGrid[-1])
        print(3)
        #if is_electron:
         #   self.ax_2.plot([energy[-1], energy[-1]], [photons[-1], 0], 'b')
          #  self.ax_2.plot([energy[-1], 10**8 + 100000], [0, 0], 'b')
        self.ax_2.set_xlim(left=-5, right=20)
        #self.ax_2.set_xlim(left=0, right= 10**8 + 100000)
        #self.ax_2.set_ylim(bottom=photons[np.where(photons==0)[0] - 1], top = max(photons)*1.05)
        self.ax_2.set_xlabel(f'Энергия фотонов, МэВ')
        print(4)
        self.ax_2.set_ylabel('Число фотонов, о.е.')
        self.ax_2.set_title(f'Спектр фотонов')
        print(5)
        self.spectrumCanvas.draw()
        print(6)

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

    def helpPopup(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("Help")
        buttonBox = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        buttonBox.clicked.connect(dlg.close)
        layout = QVBoxLayout()
        layout.addWidget(QLabel(
            "В правом меню можно переключаться между считыванием данных нейтронных датчиков и расчётом фотонейтронов"))
        layout.addWidget(QLabel(
            'В поле "Энергия электронов" указывается начальная энергия моноэнергетического пучка электронов, в правом меню выбираются единицы измерения'))
        layout.addWidget(QLabel(
            'В поле "Число электронов" указывается полное число электронов, первое число - мантисса, второе - порядок'))
        layout.addWidget(QLabel('В поле "Мишень" выбирается материал мишени'))
        layout.addWidget(QLabel('В поле "Толщина мишени" выбирается толщина мишени в метрах'))
        layout.addWidget(QLabel(
            'Одновременно с расчётом числа фотонейтронов в правом окне строится примерный вид сечения тормозного излучения'))
        layout.addWidget(buttonBox)
        dlg.setLayout(layout)
        dlg.exec()


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
        errorPopup = 1

