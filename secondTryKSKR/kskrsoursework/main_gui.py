import sys
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget,
                             QVBoxLayout, QTextEdit, QPushButton)
from PyQt6.QtGui import QFont
import pyvista as pv
from pyvistaqt import QtInteractor
import config
import meshio


# Предполагается, что у вас есть эти функции в fem_solver.py
# from fem_solver import assemble_global_stiffness, apply_bcs_and_solve

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("МКЭ Анализ Сэндвич-Панели (Ax = B)")
        self.resize(1200, 800)

        # Создаем виджет вкладок
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # Создаем сами вкладки
        self.tab_3d = QWidget()
        self.tab_math = QWidget()

        self.tabs.addTab(self.tab_3d, "Разнесенная 3D-модель (Слои)")
        self.tabs.addTab(self.tab_math, "Матрицы СЛАУ (Ax = B)")

        self.init_3d_tab()
        self.init_math_tab()

    def init_3d_tab(self):
        """Вкладка 1: Разнесенная 3D-модель слоев"""
        layout = QVBoxLayout()
        self.plotter = QtInteractor(self.tab_3d)
        layout.addWidget(self.plotter.interactor)

        btn_load = QPushButton("Загрузить и разнести слои")
        btn_load.clicked.connect(self.load_and_explode_mesh)
        layout.addWidget(btn_load)

        self.tab_3d.setLayout(layout)

    def init_math_tab(self):
        """Вкладка 2: Данные расчета Ax=B"""
        layout = QVBoxLayout()

        self.text_edit = QTextEdit()
        self.text_edit.setFont(QFont("Courier New", 10))  # Моноширинный шрифт для матриц
        self.text_edit.setReadOnly(True)
        layout.addWidget(self.text_edit)

        btn_calc = QPushButton("Вывести глобальные матрицы [A]{x}={B}")
        btn_calc.clicked.connect(self.display_matrices)
        layout.addWidget(btn_calc)

        self.tab_math.setLayout(layout)

    def load_and_explode_mesh(self):
        """Загрузка сетки и смещение слоев друг от друга для визуализации"""
        self.plotter.clear()

        # Чтение сетки (используем результаты или исходную сетку)
        mesh = pv.read(config.MESH_FILE)

        # Определяем границы слоев по оси Z (из config.py)
        z1 = config.T_OUTER
        z2 = config.T_OUTER + config.T_CORE

        # Извлекаем центры ячеек, чтобы отсортировать тетраэдры по слоям
        centers = mesh.cell_centers().points[:, 2]

        # Разделяем сетку на 3 независимых объекта
        outer_skin = mesh.extract_cells(np.where(centers <= z1)[0])
        core = mesh.extract_cells(np.where((centers > z1) & (centers <= z2))[0])
        inner_skin = mesh.extract_cells(np.where(centers > z2)[0])

        # РАЗНЕСЕНИЕ СЛОЕВ (Сдвигаем по оси Z для визуального разделения)
        shift_dist = 0.5  # Расстояние отслоения (настройте под ваши габариты)
        outer_skin.translate((0, 0, -shift_dist), inplace=True)
        inner_skin.translate((0, 0, shift_dist), inplace=True)

        # Добавляем на сцену с разными цветами
        self.plotter.add_mesh(outer_skin, show_edges=True, color="lightblue", label="Наружная обшивка")
        self.plotter.add_mesh(core, show_edges=True, color="yellow", label="Утеплитель", opacity=0.8)
        self.plotter.add_mesh(inner_skin, show_edges=True, color="lightcoral",
                              label="Внутренняя обшивка (с усилителями)")

        self.plotter.add_legend()
        self.plotter.add_axes()
        self.plotter.reset_camera()

    def display_matrices(self):
        """Формирование текстового отчета по СЛАУ: Ax = B"""
        self.text_edit.clear()
        self.text_edit.append("=== РЕШЕНИЕ СИСТЕМЫ ЛИНЕЙНЫХ АЛГЕБРАИЧЕСКИХ УРАВНЕНИЙ ===")
        self.text_edit.append("Математическая постановка: [A] {x} = {B}\n")
        self.text_edit.append(
            "Где:\n  [A] - глобальная матрица жесткости конструкции\n  {x} - вектор неизвестных перемещений узлов\n  {B} - вектор силовых нагрузок\n")

        self.text_edit.append("Выполняется сборка матриц (пожалуйста, подождите)...\n")
        QApplication.processEvents()  # Обновляем UI

        # ВНИМАНИЕ: Здесь нужно вызвать ваши реальные функции из fem_solver.py
        # K_global, points, tets = assemble_global_stiffness()
        # U, F, K_reduced = apply_bcs_and_solve(K_global, points)

        # ДЛЯ ТЕСТИРОВАНИЯ ИНТЕРФЕЙСА (заглушки случайными числами):
        # Замените эти три строки на реальные возвращаемые значения (см. выше)
        N = 15  # Количество узлов для вывода фрагмента
        K_dense_fragment = np.random.rand(N, N) * 1e8
        U_fragment = np.random.rand(N) * 1e-3
        F_fragment = np.random.rand(N) * 1000

        # ВЫВОД ФРАГМЕНТА МАТРИЦЫ [A] (K_global)
        self.text_edit.append("-" * 60)
        self.text_edit.append(f"ФРАГМЕНТ ГЛОБАЛЬНОЙ МАТРИЦЫ ЖЕСТКОСТИ [A] (Левый верхний угол {N}x{N}):")
        self.text_edit.append("-" * 60)
        for row in K_dense_fragment:
            row_str = " ".join([f"{val:10.2E}" for val in row])
            self.text_edit.append(row_str)

        # ВЫВОД ВЕКТОРА НАГРУЗОК {B} (F)
        self.text_edit.append("\n" + "-" * 60)
        self.text_edit.append(f"ВЕКТОР ВНЕШНИХ НАГРУЗОК {{B}} (Первые {N} степеней свободы):")
        self.text_edit.append("-" * 60)
        for i, val in enumerate(F_fragment):
            self.text_edit.append(f" B[{i:02d}] = {val:10.2E} Н")

        # ВЫВОД ВЕКТОРА ПЕРЕМЕЩЕНИЙ {x} (U)
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