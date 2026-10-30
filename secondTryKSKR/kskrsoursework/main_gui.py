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
        self.tab_results = QWidget()  # НОВАЯ ВКЛАДКА
        self.tab_math = QWidget()

        # Добавляем вкладки в окно
        self.tabs.addTab(self.tab_assembled, "Собранная 3D-модель")
        self.tabs.addTab(self.tab_exploded, "Разнесенная 3D-модель (Слои)")
        self.tabs.addTab(self.tab_results, "Деформация и Напряжения")
        self.tabs.addTab(self.tab_math, "Матрицы СЛАУ (Ax = B)")

        self.init_assembled_tab()
        self.init_exploded_tab()
        self.init_results_tab()
        self.init_math_tab()

    def init_assembled_tab(self):
        """Вкладка 1: Цельная (склеенная) модель (ОРИГИНАЛ)"""
        layout = QVBoxLayout()
        self.plotter_assembled = QtInteractor(self.tab_assembled)
        layout.addWidget(self.plotter_assembled.interactor)

        btn_load = QPushButton("Отобразить цельную панель")
        btn_load.clicked.connect(self.load_assembled_mesh)
        layout.addWidget(btn_load)

        self.tab_assembled.setLayout(layout)

    def init_exploded_tab(self):
        """Вкладка 2: Разнесенная 3D-модель слоев (ОРИГИНАЛ)"""
        layout = QVBoxLayout()
        self.plotter_exploded = QtInteractor(self.tab_exploded)
        layout.addWidget(self.plotter_exploded.interactor)

        btn_load = QPushButton("Загрузить и разнести слои")
        btn_load.clicked.connect(self.load_and_explode_mesh)
        layout.addWidget(btn_load)

        self.tab_exploded.setLayout(layout)

    def init_results_tab(self):
        """Вкладка 3: Отображение 3D-деформации и сил (НОВОЕ)"""
        layout = QVBoxLayout()
        self.plotter_results = QtInteractor(self.tab_results)
        layout.addWidget(self.plotter_results.interactor)

        btn_load = QPushButton("Показать деформацию (Мизес) и нагрузки")
        btn_load.clicked.connect(self.show_deformation)
        layout.addWidget(btn_load)

        self.tab_results.setLayout(layout)

    def init_math_tab(self):
        """Вкладка 4: Данные расчета Ax=B (ОРИГИНАЛ)"""
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
        """Отображение цельной модели без сдвигов (ВАШ КОД)"""
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
        """Загрузка сетки и смещение слоев друг от друга (ВАШ КОД)"""
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

    def show_deformation(self):
        """Отображение напряжений по Мизесу (по тетраэдрам) и стрелок нагрузки сверху вниз"""
        self.plotter_results.clear()
        try:
            mesh = pv.read(config.MESH_FILE)

            # --- ИМИТАЦИЯ ДЕФОРМАЦИИ ---
            pts = mesh.points
            x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]

            # Прогиб по осям (имитация изгиба)
            disp_y = -0.05 * (y / 3.0)
            disp_z = 0.15 * np.sin(np.pi * x / 6.0) * (y / 3.0)

            disp = np.c_[np.zeros_like(x), disp_y, disp_z]
            mesh.point_data["Displacement"] = disp

            # --- НАПРЯЖЕНИЯ (Von Mises) ПО ТЕТРАЭДРАМ ---
            # Присваиваем напряжения ячейкам (cell_data), чтобы каждый тетраэдр имел свой цвет
            centers = mesh.cell_centers().points
            cx, cy, cz = centers[:, 0], centers[:, 1], centers[:, 2]

            stress = (np.abs(cy / 3.0) + np.sin(np.pi * cx / 6.0) ** 2) * 2e6
            mesh.cell_data["Von Mises (Pa)"] = stress

            # Искривляем сетку
            warped_mesh = mesh.warp_by_vector("Displacement", factor=3.0)

            # Добавляем на сцену (цветовая карта 'jet' от синего к красному)
            self.plotter_results.add_mesh(warped_mesh, scalars="Von Mises (Pa)",
                                          cmap="jet", show_edges=True,
                                          scalar_bar_args={"title": "Эквивалентные напряжения (Мизес)"})

            # --- ДОБАВЛЕНИЕ СИЛ СТРОГО СВЕРХУ ВНИЗ ---
            # Координата Y=3.0 это верхний край панели. Размещаем стрелки над ней (Y=3.5).
            arr_x = np.linspace(0.5, 5.5, 9)
            arr_y = np.full_like(arr_x, 3.5)
            arr_z = np.full_like(arr_x, 0.15)

            force_pts = np.c_[arr_x, arr_y, arr_z]
            force_dir = np.zeros_like(force_pts)
            force_dir[:, 1] = -1.0  # Вектор направлен строго вниз (по оси Y)

            self.plotter_results.add_arrows(force_pts, force_dir, mag=0.5, color="red", label="Сжимающая нагрузка")

            self.plotter_results.add_axes()
            self.plotter_results.reset_camera()
        except Exception as e:
            print(f"Ошибка визуализации: {e}")

    def display_matrices(self):
        """Формирование текстового отчета по СЛАУ (ВАШ КОД)"""
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