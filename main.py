from GUI import *

if __name__ == '__main__':
    app = QApplication(sys.argv)
    photon = PhotonWindow()
    photon.show()
    sys.exit(app.exec())

