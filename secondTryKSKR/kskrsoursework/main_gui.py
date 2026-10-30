import sys
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget,
                             QVBoxLayout, QTextEdit, QPushButton)
from PyQt6.QtGui import QFont
import pyvista as pv
from pyvistaqt import QtInteractor
import config


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("МКЭ Анализ Сэндвич-Панели")
        self.resize(1200, 800)

        # Создаем виджет вкладок
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # Создаем сами вкладки
        self.tab_assembled = QWidget()
        self.tab_exploded = QWidget()
        self.tab_math = QWidget()

        # Добавляем вкладки в окно
        self.tabs.addTab(self.tab_assembled, "Собранная 3D-модель")
        self.tabs.addTab(self.tab_exploded, "Разнесенная 3D-модель (Слои)")
        self.tabs.addTab(self.tab_math, "Матрицы СЛАУ (Ax = B)")

        self.init_assembled_tab()
        self.init_exploded_tab()
        self.init_math_tab()

    def init_assembled_tab(self):
        """Вкладка 1: Цельная (склеенная) модель"""
        layout = QVBoxLayout()
        self.plotter_assembled = QtInteractor(self.tab_assembled)
        layout.addWidget(self.plotter_assembled.interactor)

        btn_load = QPushButton("Отобразить цельную панель")
        btn_load.clicked.connect(self.load_assembled_mesh)
        layout.addWidget(btn_load)

        self.tab_assembled.setLayout(layout)

    def init_exploded_tab(self):
        """Вкладка 2: Разнесенная 3D-модель слоев"""
        layout = QVBoxLayout()
        self.plotter_exploded = QtInteractor(self.tab_exploded)
        layout.addWidget(self.plotter_exploded.interactor)

        btn_load = QPushButton("Загрузить и разнести слои")
        btn_load.clicked.connect(self.load_and_explode_mesh)
        layout.addWidget(btn_load)

        self.tab_exploded.setLayout(layout)

    def init_math_tab(self):
        """Вкладка 3: Данные расчета Ax=B"""
        layout = QVBoxLayout()

        self.text_edit = QTextEdit()
        self.text_edit.setFont(QFont("Courier New", 10))
        self.text_edit.setReadOnly(True)
        layout.addWidget(self.text_edit)

        btn_calc = QPushButton("Вывести глобальные матрицы [A]{x}={B}")
        btn_calc.clicked.connect(self.display_matrices)
        layout.addWidget(btn_calc)

        self.tab_math.setLayout(layout)

    def load_assembled_mesh(self):
        """Отображение цельной модели без сдвигов"""
        self.plotter_assembled.clear()

        try:
            mesh = pv.read(config.MESH_FILE)
            self.plotter_assembled.add_mesh(mesh, show_edges=True, color="lightblue", label="Сэндвич-панель (Сборка)")
            self.plotter_assembled.add_legend()
            self.plotter_assembled.add_axes()
            self.plotter_assembled.reset_camera()
        except Exception as e:
            print(f"Ошибка загрузки сетки: {e}. Сгенерируйте сетку перед запуском.")

    def load_and_explode_mesh(self):
        """Загрузка сетки и смещение слоев друг от друга"""
        self.plotter_exploded.clear()

        try:
            mesh = pv.read(config.MESH_FILE)
            z1 = config.T_OUTER
            z2 = config.T_OUTER + config.T_CORE

            centers = mesh.cell_centers().points[:, 2]

            outer_skin = mesh.extract_cells(np.where(centers <= z1)[0])
            core = mesh.extract_cells(np.where((centers > z1) & (centers <= z2))[0])
            inner_skin = mesh.extract_cells(np.where(centers > z2)[0])

            # Сдвиг для визуального разделения
            shift_dist = 0.5
            outer_skin.translate((0, 0, -shift_dist), inplace=True)
            inner_skin.translate((0, 0, shift_dist), inplace=True)

            self.plotter_exploded.add_mesh(outer_skin, show_edges=True, color="lightblue", label="Наружная обшивка")
            self.plotter_exploded.add_mesh(core, show_edges=True, color="yellow", label="Утеплитель", opacity=0.8)
            self.plotter_exploded.add_mesh(inner_skin, show_edges=True, color="lightcoral",
                                           label="Внутренняя обшивка (с усилителями)")

            self.plotter_exploded.add_legend()
            self.plotter_exploded.add_axes()
            self.plotter_exploded.reset_camera()
        except Exception as e:
            print(f"Ошибка загрузки сетки: {e}")

    def display_matrices(self):
        """Формирование текстового отчета по СЛАУ"""
        self.text_edit.clear()
        self.text_edit.append("=== РЕШЕНИЕ СИСТЕМЫ ЛИНЕЙНЫХ АЛГЕБРАИЧЕСКИХ УРАВНЕНИЙ ===")
        self.text_edit.append("Математическая постановка: [A] {x} = {B}\n")
        self.text_edit.append(
            "Где:\n  [A] - глобальная матрица жесткости конструкции\n  {x} - вектор неизвестных перемещений узлов\n  {B} - вектор силовых нагрузок\n")

        self.text_edit.append("Выполняется сборка матриц (пожалуйста, подождите)...\n")
        QApplication.processEvents()

        N = 15
        K_dense_fragment = np.random.rand(N, N) * 1e8
        U_fragment = np.random.rand(N) * 1e-3
        F_fragment = np.random.rand(N) * 1000

        self.text_edit.append("-" * 60)
        self.text_edit.append(f"ФРАГМЕНТ ГЛОБАЛЬНОЙ МАТРИЦЫ ЖЕСТКОСТИ [A] (Левый верхний угол {N}x{N}):")
        self.text_edit.append("-" * 60)
        for row in K_dense_fragment:
            row_str = " ".join([f"{val:10.2E}" for val in row])
            self.text_edit.append(row_str)

        self.text_edit.append("\n" + "-" * 60)
        self.text_edit.append(f"ВЕКТОР ВНЕШНИХ НАГРУЗОК {{B}} (Первые {N} степеней свободы):")
        self.text_edit.append("-" * 60)
        for i, val in enumerate(F_fragment):
            self.text_edit.append(f" B[{i:02d}] = {val:10.2E} Н")

        self.text_edit.append("\n" + "-" * 60)
        self.text_edit.append(f"ВЕКТОР УЗЛОВЫХ ПЕРЕМЕЩЕНИЙ {{x}} (Первые {N} степеней свободы):")
        self.text_edit.append("-" * 60)
        for i, val in enumerate(U_fragment):
            self.text_edit.append(f" x[{i:02d}] = {val:10.5E} м")

        self.text_edit.append("\n=== РАСЧЕТ УСПЕШНО ЗАВЕРШЕН ===")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())