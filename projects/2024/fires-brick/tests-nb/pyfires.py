#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Engy-4390 course
import numpy as np 
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def get_domain_partition(degree, n_elem, x_min, x_max, bc_x_min='essential', bc_x_max='essential'):
    #assert degree == 1
    # Local node numbering on parent domain
    # --0--------------1---->
    #  -1      0      +1    zetta
    gnodes_x = np.linspace(x_min, x_max, n_elem*degree+1, dtype=np.float64)
    patches = list()
    local_to_global_node_id_map = list()
    for e in range(n_elem):
        gnode_id_1 = degree*e   #left
        gnode_id_2 = degree*e+1 #center
        gnode_id_3 = degree*e+2 #right
        x1 = gnodes_x[gnode_id_1]
        x2 = gnodes_x[gnode_id_2]
        x3 = gnodes_x[gnode_id_3]
        # Local node id:  0   1
        patches.append((x1, x2, x3))
        # Local node id:                        0           2             1
        local_to_global_node_id_map.append([gnode_id_1, gnode_id_2,  gnode_id_3])
    if bc_x_min == 'essential':
        local_to_global_node_id_map[0][0] = -1
    if bc_x_max == 'essential':
        local_to_global_node_id_map[-1][-1] = -1
    return (patches, gnodes_x, local_to_global_node_id_map)

def get_parent_mapping():
    # zetta in [-1,1]
    parent_mapping = lambda zetta, x_e_bar, h_e: x_e_bar + h_e/2*zetta # compute x
    parent_mapping_prime = lambda h_e: h_e/2                           # compute mapping derivative wrt zetta
    # x in Omega_e
    inverse_parent_mapping = lambda x, x_e_bar, h_e: (x - x_e_bar)*2/h_e # compute zetta
    return (parent_mapping, parent_mapping_prime, inverse_parent_mapping)

def get_parent_basis_functions():
    parent_basis_func_list = list()
    parent_basis_func_prime_list = list()
    parent_basis_func_list.append(lambda zetta: (zetta**2-zetta)/2)  # left
    parent_basis_func_list.append(lambda zetta: -(zetta**2-1))  # middle
    parent_basis_func_list.append(lambda zetta:  (zetta**2+zetta)/2)  # right
    parent_basis_func_prime_list.append(lambda zetta: (2*zetta-1)/2) # left
    parent_basis_func_prime_list.append(lambda zetta: -2*zetta) # middle
    parent_basis_func_prime_list.append(lambda zetta:  (2*zetta+1)/2) # right
    return (parent_basis_func_list, parent_basis_func_prime_list)

def global_basis_function(i, x, domain_partition, parent_mapping, parent_basis_functions):
    try:
        len(x)
    except TypeError:
        x = np.array([x])
    if not isinstance(x, np.ndarray):
       assert isinstance(x, list) or isinstance(x, tuple)
       x = np.array(x)
    phi_i_x = np.copy(x) * 0.0 # initialization
    phi_prime_i_x = np.copy(x) * 0.0 # initialization
    patches = domain_partition[0]
    local_to_global_node_id_map = domain_partition[2]
    inverse_parent_mapping = parent_mapping[2]
    parent_basis_func_list = parent_basis_functions[0]
    parent_basis_func_prime_list = parent_basis_functions[1]
    # expensive reverse lookup
    for j, x_j in enumerate(x):
        for e, nodes_x in enumerate(patches):
            if nodes_x[0] <= x_j <= nodes_x[-1]:
                n_lnodes = len(nodes_x)
                for I in range(n_lnodes):
                    if local_to_global_node_id_map[e][I] == i:
                        x_e_bar = (nodes_x[0] + nodes_x[-1])/2
                        h_e = nodes_x[-1] - nodes_x[0]
                        zetta = inverse_parent_mapping(x_j, x_e_bar, h_e)
                        phi_i_x[j] = parent_basis_func_list[I](zetta)
                        phi_prime_i_x[j] = parent_basis_func_prime_list[I](zetta)
                break
    return [phi_i_x, phi_prime_i_x]

def get_global_basis_functions(domain_partition, parent_mapping, parent_basis_functions,global_basis_function):
    basis_func_list = list()
    basis_func_prime_list = list()
    n_elem = len(domain_partition[0])
    n_gnodes = domain_partition[1].size
    local_to_global_node_id_map = domain_partition[2]
    phi_i = lambda i, x: global_basis_function(i, x, domain_partition, parent_mapping, parent_basis_functions)[0]
    phi_prime_i = lambda i, x: global_basis_function(i, x, domain_partition,parent_mapping, parent_basis_functions)[1]
    visited = [False]*n_gnodes
    for e in range(n_elem):
        for I in range(len(local_to_global_node_id_map[e])):
            gnode_id = local_to_global_node_id_map[e][I]
            if gnode_id >= 0 and not visited[gnode_id]:
                      basis_func_list.append(lambda x, i=gnode_id: phi_i(i,x))
                      basis_func_prime_list.append(lambda x, i=gnode_id: phi_prime_i(i,x))

                      visited[gnode_id] = True
    assert len(basis_func_list) >= 1, 'There are no basis functions to build.'
    return [basis_func_list, basis_func_prime_list]

def inner_product(u, v, patches):
    integrand = lambda x: u(x) * v(x)
    inner_product = 0.0
    epsrel = 1e-6
    epsabs = 1e-6
    limit = 150
    for nodes_x in patches:
        (inner_product_e, _) = quad(integrand, nodes_x[0], nodes_x[-1],
                                    epsrel=epsrel, epsabs=epsabs, limit=limit)
        inner_product += inner_product_e
    return inner_product

def linear_func(x_min, x_max, vals):
    """1D function with interpolation wrapper.
    """

    cond_shape_pts = [(x_min, vals[0]), (x_max, vals[1])]
    cond_shape_pts_mtrx = np.array(cond_shape_pts)
    cte_func = interp1d(cond_shape_pts_mtrx[:,0], cond_shape_pts_mtrx[:,1])

    return cte_func

def plot_func(func, x_min, x_max, n_pts, xlabel='xlabel', ylabel='ylabel', title='title',
              x_scale=1.0, y_scale=1.0, gold_data=None):

    if gold_data is not None:
        assert isinstance(gold_data, np.ndarray) and len(gold_data.shape) == 2

    plt.figure()
    x_vals = np.linspace(x_min, x_max, n_pts)
    if gold_data is not None:
        plt.plot(x_vals*x_scale, func(x_vals)*y_scale, '-b', gold_data[:,0], gold_data[:,1], '--r')
    else:
        plt.plot(x_vals*x_scale, func(x_vals)*y_scale)
    plt.title (title)
    plt.xlabel(xlabel)
    plt.ylabel (ylabel)
    plt.grid()
    plt.show()

def u_star(x, phi_lst, lift_func, c_vec):
    g_x = lift_func(x)
    for (j, phi_i) in enumerate(phi_lst):
        g_x = g_x + c_vec[j] * phi_i(x)
    return g_x

def u_star_prime(x, phi_prime_lst, lift_func_prime, c_vec, h_e):
    g_x = lift_func_prime(x)
    for j in range(len(phi_prime_lst)):
        g_x = g_x + c_vec[j] * 2/h_e * phi_prime_lst[j](x)
    return g_x

def build_a_mtrx(phi_lst, phi_prime_lst, k_func, domain_partition, x_min, x_max, n_elem):

    A_mtrx = np.zeros((len(phi_lst), len(phi_lst)), dtype=np.float64)
    patches = domain_partition[0]
    for i in range(len(phi_lst)):
        for j in range(len(phi_lst)):

            phi_i=phi_lst[i]
            phi_j=phi_lst[j]

            phi_prime_i=phi_prime_lst[i]
            phi_prime_j=phi_prime_lst[j]

            h_e=(x_max-x_min)/n_elem

            d_x_phi_prime_j = lambda x: k_func(x) * ((2/h_e)*phi_prime_j(x))

            prima = lambda x: phi_prime_i(x)*(2/h_e)

            A_mtrx[i,j] = inner_product(prima, d_x_phi_prime_j, patches)

    return A_mtrx

def build_b_vec(phi_list, phi_prime_list,
                k_func, f_func, lift_func_prime, domain_partition, x_min, x_max, n_elem):

    b_vec = np.zeros(len(phi_list), dtype=np.float64)
    patches = domain_partition[0]

    for i in range(len(phi_list)):
        phi_i=phi_list[i]
        phi_prime_i=phi_prime_list[i]

        h_e=(x_max-x_min)/n_elem

        b_vec[i] = inner_product(f_func, phi_i, patches)

        first_term = lambda x: lift_func_prime(x) * k_func(x)
        phi_prima_i = lambda x: phi_prime_i(x)*(2/h_e)

        b_vec[i] -= inner_product(first_term, phi_prima_i, patches)

    return b_vec

def build_b_vec_2(phi_list, phi_prime_list, 
                  k_func, f_func, lift_func_prime, domain_partition, x_min, x_max, n_elem, htc, u_a, u_b):
    
    b_vec_2 = np.zeros(len(phi_list), dtype=np.float64)
    patches = domain_partition[0]

    for i in range(len(phi_list)):
        phi_i=phi_list[i]
        phi_prime_i=phi_prime_list[i]

        h_e=(x_max-x_min)/n_elem
    
        b_vec_2[i] = inner_product(f_func, phi_i, patches)-htc*(f_func(x_max)-u_b)*phi_i(x_max)-htc*(f_func(x_min)-u_a)*phi_i(x_min)
    
        first_term = lambda x: lift_func_prime(x)*k_func(x)
        phi_prima_i = lambda x: phi_prime_i(x)*(2/h_e)
    
        b_vec_2[i] -= inner_product(first_term, phi_prima_i, patches)
        
    return b_vec_2

    