# -*- coding:utf-8 -*-

import numpy as np
import pickle
import itertools
import random

mol_interact_relations = pickle.load(open('../Data/data1_molecular_link_index.pkl', 'rb'))
mol_index = pickle.load(open('../Data/data2_smiles_index.pkl', 'rb'))
mol_syno = pickle.load(open('../Data/data3_synonym_smiles.pkl', 'rb'))
can_mol_index = pickle.load(open('../Data/data4_smiles_index.pkl', 'rb'))
mol_fps = pickle.load(open('../Data/data5_fingerprint_3d.pkl', 'rb'))

substr_weight = 0.02
molecule_num = 2947
can_mol_dict = dict(can_mol_index)
input_smi_index = can_mol_dict
mol_input_fps = mol_fps

# 建立网络的边并赋予不同的权重
substructure_matrix = np.array(mol_input_fps, dtype=np.float64)
substructure_matrix = substructure_matrix[:, np.sum(substructure_matrix, axis=0) != 0]  # 删除全为0的特征列
mol_num, substructure_num = substructure_matrix.shape  # 得到行数、列数
substructure_links = []
for mol in range(mol_num):
    for i in range(substructure_num):
        if substructure_matrix[mol, i] == 1:
            substructure_links.append([mol, mol_num + i])
substructure_links = [item + [substr_weight] for item in substructure_links]
mol_interact_relations = [item + [1 / substr_weight] for item in mol_interact_relations]
links = mol_interact_relations + substructure_links

# 网络矩阵的建立
mat_nodes = list(itertools.chain.from_iterable(links))
mat_nodes = set(mat_nodes)
mat_nodes = {np.int32(node) for node in mat_nodes}
mat_size = np.int32(max(mat_nodes) + 1)
network = np.zeros((mat_size, mat_size))
for item in links:
    network[np.int32(item[0]), np.int32(item[1])] = item[2]
b_network = network.copy()   #半网络
network = network + network.T

# 将网络连接矩阵转换为得分矩阵
degree = np.tile(np.sum(network, 0), (np.size(network, 0), 1))
degree = degree.T
trans_mat = np.divide(network, degree)
unit_mat = np.eye(np.size(trans_mat, 0))
score = unit_mat
step_i = 0
while step_i < 3:
    score = np.dot(trans_mat.T, score)
    step_i += 1
score = score + score.T
score[np.isnan(score)] = 0
score[np.isinf(score)] = 0

####后续为网络评价相关内容
# 提取0-2946范围内的节点构成的子矩阵
sub_matrix = network[:molecule_num, :molecule_num]
b_sub_matrix = b_network[:molecule_num, :molecule_num]
# 找到network矩阵中值为0的边的位置
zero_edges = np.argwhere(sub_matrix == 0)    # 不存在边的位置
# 现在改为之后进行收集

#获取存在的链接的索引
b_real_edges = np.argwhere(b_sub_matrix == 1/substr_weight)  # 得到real edges坐标,注意，因为网络是对称的，所以只需要获取b_sub_matrix中的


f_average_results = []
for iteration in range(10):   # 执行十次十折
    num_folds = 10
    num_edges = len(b_real_edges)
    results = []

    np.random.shuffle(b_real_edges)  # 打乱真实边的坐标数据
    # 将数据划分成十折
    fold_size = len(b_real_edges) // 10
    folds = [b_real_edges[i:i + fold_size] for i in range(0, len(b_real_edges), fold_size)]  # 包含十折数据

    # 执行十折交叉验证
    for fold in range(num_folds):
        # 从所有索引中选择一份
        validation_data = folds[fold]
        random_edges = validation_data
        # 将这部分坐标的权重值设置为0，然后将b_network复制，即可得到对称的去掉相应位置的整个网络
        b_network_copy = b_network.copy()
        for edge in random_edges:
            b_network_copy[tuple(edge)] = 0
        new_network = b_network_copy + b_network_copy.T   # 此网络为去除随机边后的新网络

        # 将新网络连接矩阵转换为得分矩阵
        degree = np.tile(np.sum(new_network, 0), (np.size(new_network, 0), 1))
        degree = degree.T
        trans_mat = np.divide(new_network, degree)
        unit_mat = np.eye(np.size(trans_mat, 0))
        score = unit_mat
        step_i = 0
        while step_i < 3:
            score = np.dot(trans_mat.T, score)
            step_i += 1
        score = score + score.T
        score[np.isnan(score)] = 0
        score[np.isinf(score)] = 0

        # 收集对应位置在score矩阵中的值
        test_edges_values = [score[row, col] for row, col in random_edges]
        test_edges_list = list(test_edges_values)
        null_edges_values = [score[row, col] for row, col in zero_edges]
        null_edges_list = list(null_edges_values)

        count_A = 0
        count_B = 0
        for _ in range(10000):
            num1 = random.choice(test_edges_list)
            num2 = random.choice(null_edges_list)
            if num1 > num2:
                count_A += 1
            elif num1 == num2:
                count_B += 1
        result = (count_A + 0.5 * count_B) / 10000
        print(f"Fold {fold + 1} Result: {result}")
        results.append(result)
    average_result = np.mean(results)
    print(f"Average Result: {average_result}")
    print('完成一次十折交叉验证')
    print('----------')
    f_average_results.append(average_result)

final_average_result = np.mean(f_average_results)
print(f"\nFinal Average Result: {final_average_result}")

# 计算十次十折交叉验证的标准差
final_std_dev = np.std(f_average_results)
print(f"\nFinal Standard Deviation: {final_std_dev}")