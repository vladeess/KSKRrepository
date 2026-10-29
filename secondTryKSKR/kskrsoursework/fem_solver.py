import numpy as np
import meshio
import scipy.sparse as sp
import config
import scipy.sparse.linalg as spla

import numpy as np
import meshio
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import config


def print_formatted_matrix(name, matrix):
    """Красивый вывод матриц в консоль для отчета"""
    print(f"\n--- {name} ({matrix.shape[0]}x{matrix.shape[1]}) ---")
    for row in matrix:
        print("  " + " ".join([f"{val:11.3E}" for val in row]))


# ... (функции get_elasticity_matrix и calculate_tet_volume_and_B оставляем без изменений)

def assemble_global_stiffness():
    print("\n[1/4] Чтение сетки из файла...")
    mesh = meshio.read(config.MESH_FILE)
    points = mesh.points
    num_nodes = len(points)
    total_dofs = num_nodes * 3
    print(f"  Узлов: {num_nodes}, Степеней свободы: {total_dofs}")

    row_indices, col_indices, values = [], [], []

    D_outer = get_elasticity_matrix(config.MATERIALS["outer_skin"]["E"], config.MATERIALS["outer_skin"]["nu"])
    D_core = get_elasticity_matrix(config.MATERIALS["core"]["E"], config.MATERIALS["core"]["nu"])
    D_inner = get_elasticity_matrix(config.MATERIALS["inner_skin"]["E"], config.MATERIALS["inner_skin"]["nu"])

    # ВЫВОД МАТРИЦ УПРУГОСТИ ДЛЯ ОТЧЕТА
    print("\n================ РАСЧЕТЫ ДЛЯ ПОЯСНИТЕЛЬНОЙ ЗАПИСКИ ================")
    print_formatted_matrix("Матрица упругости [D] (Обшивка - Сталь)", D_outer)
    print_formatted_matrix("Матрица упругости [D] (Утеплитель)", D_core)

    tets_list = []
    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            tets_list.append(cell_block.data)
    if not tets_list:
        raise ValueError("В сетке не найдены тетраэдры!")
    tets = np.vstack(tets_list)

    print(f"\n[3/4] Вычисление локальных матриц и ансамблирование ({len(tets)} элементов)...")

    for count, tet in enumerate(tets):
        tet_nodes = points[tet]

        z_center = np.mean(tet_nodes[:, 2])
        if z_center <= config.T_OUTER:
            D = D_outer
        elif z_center <= (config.T_OUTER + config.T_CORE):
            D = D_core
        else:
            D = D_inner

        vol, B = calculate_tet_volume_and_B(tet_nodes)
        K_e = B.T @ D @ B * vol

        # ВЫВОД ПОДРОБНОГО РАСЧЕТА ДЛЯ ПЕРВОГО ЭЛЕМЕНТА
        if count == 0:
            print(f"\n--- Детальный расчет для Элемента №0 (Узлы: {tet}) ---")
            print(f"  Объем тетраэдра V = {vol:.6E} м^3")
            print_formatted_matrix("Градиентная матрица деформаций [B]", B)
            print_formatted_matrix("Локальная матрица жесткости [K_e]", K_e)
            print("===================================================================\n")

        dofs = np.zeros(12, dtype=int)
        for i in range(4):
            dofs[i * 3] = tet[i] * 3
            dofs[i * 3 + 1] = tet[i] * 3 + 1
            dofs[i * 3 + 2] = tet[i] * 3 + 2

        for i in range(12):
            for j in range(12):
                row_indices.append(dofs[i])
                col_indices.append(dofs[j])
                values.append(K_e[i, j])

    print("[4/4] Создание глобальной разреженной матрицы...")
    K_global = sp.coo_matrix((values, (row_indices, col_indices)), shape=(total_dofs, total_dofs)).tocsr()
    return K_global, points, tets


# ... (Остальная часть файла с apply_bcs_and_solve остается без изменений)

def get_elasticity_matrix(E, nu):
    """Матрица упругости [D] (6x6) для изотропного материала в 3D"""
    factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    D = np.zeros((6, 6))
    D[0, 0] = D[1, 1] = D[2, 2] = 1.0 - nu
    D[0, 1] = D[1, 0] = D[0, 2] = D[2, 0] = D[1, 2] = D[2, 1] = nu
    D[3, 3] = D[4, 4] = D[5, 5] = (1.0 - 2.0 * nu) / 2.0
    return factor * D


def calculate_tet_volume_and_B(nodes):
    """Вычисление объема тетраэдра и матрицы [B] (6x12)"""
    # Матрица координат узлов для вычисления определителя (Якобиана)
    ones = np.ones((4, 1))
    C = np.hstack((ones, nodes))

    volume = np.linalg.det(C) / 6.0

    # Если порядок обхода узлов дает отрицательный объем, берем модуль
    # (в идеале Gmsh выдает правильный порядок, но перестраховка не помешает)
    abs_volume = abs(volume)

    ic = np.linalg.inv(C)

    # Производные функций формы по x, y, z (это строки инвертированной матрицы C, начиная со 2-й)
    dN_dx = ic[1, :]
    dN_dy = ic[2, :]
    dN_dz = ic[3, :]

    B = np.zeros((6, 12))
    for i in range(4):
        B[0, i * 3] = dN_dx[i]
        B[1, i * 3 + 1] = dN_dy[i]
        B[2, i * 3 + 2] = dN_dz[i]

        B[3, i * 3] = dN_dy[i]
        B[3, i * 3 + 1] = dN_dx[i]

        B[4, i * 3 + 1] = dN_dz[i]
        B[4, i * 3 + 2] = dN_dy[i]

        B[5, i * 3] = dN_dz[i]
        B[5, i * 3 + 2] = dN_dx[i]

    return abs_volume, B


def assemble_global_stiffness():
    print("[1/4] Чтение сетки из файла...")
    mesh = meshio.read(config.MESH_FILE)
    points = mesh.points
    num_nodes = len(points)
    total_dofs = num_nodes * 3

    print(f"  Узлов: {num_nodes}, Степеней свободы: {total_dofs}")

    # Списки для сборки разреженной матрицы (формат COO)
    row_indices = []
    col_indices = []
    values = []

    print("[2/4] Подготовка свойств материалов...")
    # Маппинг материалов (для упрощения примера предполагаем,
    # что тетраэдры идут в том же порядке, что и физические группы Gmsh)
    # В реальной задаче нужно читать cell_data["gmsh:physical"], если Gmsh экспортировал их.
    # Здесь мы определяем слой по Z-координате центроида элемента.

    D_outer = get_elasticity_matrix(config.MATERIALS["outer_skin"]["E"], config.MATERIALS["outer_skin"]["nu"])
    D_core = get_elasticity_matrix(config.MATERIALS["core"]["E"], config.MATERIALS["core"]["nu"])
    D_inner = get_elasticity_matrix(config.MATERIALS["inner_skin"]["E"], config.MATERIALS["inner_skin"]["nu"])

    # ... (после строк с D_inner = ...)
    tets_list = []
    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            tets_list.append(cell_block.data)

    if not tets_list:
        raise ValueError("В сетке не найдены тетраэдры!")

    # Объединяем тетраэдры всех слоев в одну матрицу
    tets = np.vstack(tets_list)

    print(f"[3/4] Вычисление локальных матриц и ансамблирование ({len(tets)} элементов)...")
    # ... (дальше идет цикл for count, tet in enumerate(tets):)

    for count, tet in enumerate(tets):
        tet_nodes = points[tet]

        # Определяем материал по координате Z центроида
        z_center = np.mean(tet_nodes[:, 2])
        if z_center <= config.T_OUTER:
            D = D_outer
        elif z_center <= (config.T_OUTER + config.T_CORE):
            D = D_core
        else:
            D = D_inner

        vol, B = calculate_tet_volume_and_B(tet_nodes)

        # Локальная матрица жесткости: Ke = B^T * D * B * V
        K_e = B.T @ D @ B * vol

        # Глобальные индексы степеней свободы для 4 узлов элемента
        dofs = np.zeros(12, dtype=int)
        for i in range(4):
            dofs[i * 3] = tet[i] * 3
            dofs[i * 3 + 1] = tet[i] * 3 + 1
            dofs[i * 3 + 2] = tet[i] * 3 + 2

        # Заполнение векторов для COO матрицы
        for i in range(12):
            for j in range(12):
                row_indices.append(dofs[i])
                col_indices.append(dofs[j])
                values.append(K_e[i, j])

        if count % 5000 == 0 and count > 0:
            print(f"  Обработано {count} элементов...")

    print("[4/4] Создание глобальной разреженной матрицы...")
    K_global = sp.coo_matrix((values, (row_indices, col_indices)), shape=(total_dofs, total_dofs)).tocsr()

    return K_global, points, tets





def apply_bcs_and_solve(K_global, points):
    print("[1/3] Поиск граничных узлов и приложения нагрузки...")
    num_nodes = len(points)
    total_dofs = num_nodes * 3

    # Индексы степеней свободы для заделки и узлы для нагрузки
    fixed_dofs = []
    loaded_nodes = []

    TOL = 1e-4  # Допуск для поиска координат с плавающей точкой
    for i, (x, y, z) in enumerate(points):
        # Жесткая заделка снизу (грань Y = 0)
        if abs(y) < TOL:
            # Блокируем перемещения u, v, w (все 3 степени свободы узла)
            fixed_dofs.extend([i * 3, i * 3 + 1, i * 3 + 2])

        # Поиск узлов на верхней грани для приложения нагрузки (Y = B)
        if abs(y - config.B) < TOL:
            loaded_nodes.append(i)

    print(f"  Зафиксировано степеней свободы: {len(fixed_dofs)}")
    print(f"  Узлов под нагрузкой: {len(loaded_nodes)}")

    # Формируем глобальный вектор нагрузок
    F_full = np.zeros(total_dofs)
    total_force_y = -50000.0  # Суммарная нагрузка в Ньютонах (-50 кН вниз)

    if len(loaded_nodes) > 0:
        # Упрощенное равномерное распределение силы по верхним узлам
        force_per_node = total_force_y / len(loaded_nodes)
        for node in loaded_nodes:
            F_full[node * 3 + 1] = force_per_node  # +1 соответствует оси Y

    print("[2/3] Урезание матриц с учетом граничных условий (слайсинг)...")
    all_dofs = np.arange(total_dofs)
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    # Вырезаем только свободные степени свободы из глобальной матрицы и вектора сил
    K_reduced = K_global[free_dofs, :][:, free_dofs]
    F_reduced = F_full[free_dofs]

    print("[3/3] Решение СЛАУ...")
    # Решение системы K*U = F для свободных узлов
    U_reduced = spla.spsolve(K_reduced, F_reduced)

    # Восстановление полного вектора перемещений (включая нули в заделке)
    U_full = np.zeros(total_dofs)
    U_full[free_dofs] = U_reduced

    max_disp = np.max(np.abs(U_full)) * 1000  # перевод в мм
    print("  Решение успешно найдено!")
    print(f"  Максимальное перемещение: {max_disp:.3f} мм")

    return U_full


# --- Обновите блок запуска в самом низу файла ---
if __name__ == "__main__":
    K_global, points, tets = assemble_global_stiffness()
    print("Матрица жесткости собрана. Размер:", K_global.shape)

    U_full = apply_bcs_and_solve(K_global, points)


if __name__ == "__main__":
    K, pts, tets = assemble_global_stiffness()
    print("Матрица успешно собрана! Размер:", K.shape)