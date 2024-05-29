from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QFileDialog, QVBoxLayout, QWidget
import sys

class MainWindow(QMainWindow):
    def init(self, parent=None):
        super(MainWindow, self).init(parent)
        self.setWindowTitle('2D FFT GUI')

        # Создаем кнопку для загрузки файла
        self.load_file_button = QPushButton('Load File')
        self.load_file_button.clicked.connect(self.load_file)

        # Создаем место для отображения имени загруженного файла
        self.file_label = QLabel()

        # Размещаем кнопку и метку вертикально
        layout = QVBoxLayout()
        layout.addWidget(self.load_file_button)
        layout.addWidget(self.file_label)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def load_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_name, _ = QFileDialog.getOpenFileName(self,'QFileDialog.getOpenFileName()', '','All Files (*);;Python Files (*.py)', options=options)
        if file_name:
            self.file_label.setText(file_name)

if name == 'main':
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())
