from GUI import *

def onChanged():
    if w.currentIndex() == 0:
        w.setCurrentIndex(progs[counter.ProgramChannelBox.currentText()])
    elif w.currentIndex() == 1:
        w.setCurrentIndex(progs[photon.ProgramChannelBox.currentText()])
    elif w.currentIndex() == 2:
        w.setCurrentIndex(progs[journal.ProgramChannelBox.currentText()])

if __name__ == '__main__':
    app = QApplication(sys.argv)
    counter = CounterWindow()
    photon = PhotonWindow()
    journal = JournalWindow()
    w = QStackedWidget()
    w.addWidget(counter)

    w.addWidget(photon)
    w.addWidget(journal)
    w.resize(1400,950)
    #w.setWindowTitle('MultiNeutronator')
    counter.ProgramChannelBox.activated.connect(onChanged)
    photon.ProgramChannelBox.activated.connect(onChanged)
    journal.ProgramChannelBox.activated.connect(onChanged)
    w.show()
    sys.exit(app.exec())

