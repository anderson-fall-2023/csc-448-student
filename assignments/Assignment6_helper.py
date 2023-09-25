import numpy as np
import pandas as pd

def greedy_sorting(P):
    P = P.copy()
    approx_rev_distance = 0
    ## BEGIN SOLUTION
    for k in P.index:
        if P.loc[k] != k:
            loc = int(np.where((P == k) | (P == -k))[0])
            P.loc[k:(loc+1)] = list(-1*P.loc[k:(loc+1)])[::-1]
                
            approx_rev_distance += 1
            if P.loc[k] == -k:
                P.loc[k] = k
                approx_rev_distance += 1
    ## END SOLUTION
    # YOUR SOLUTION HERE
    return approx_rev_distance

def count_breakpoints(P):
    nbreakpoints = 0
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    if P.iloc[0] != P.index[0]:
        nbreakpoints += 1
    if P.iloc[-1] != P.index[-1]:
        nbreakpoints += 1    
    for i in range(len(P)-1):
        if P.iloc[i+1] - P.iloc[i] != 1:
            nbreakpoints += 1
    ## END SOLUTION
    return nbreakpoints

import networkx as nx
import copy

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

def to_adj(T):
    df = pd.DataFrame(nx.adjacency_matrix(T).todense(),index=T.nodes(),columns=T.nodes())
    return df

def to_edge_list(T):
    return list(T.edges())

def show(G,G2=None,P_G=None,red_blue=True):
    if P_G is None:
        pos = nx.circular_layout(G)
    else:
        pos = nx.circular_layout(P_G)            

    nx.draw_networkx_nodes(G, pos=pos,node_size=600,node_color='w')
    nx.draw_networkx_labels(G, pos=pos)        
    nodes = list(G.nodes())
    edge_list_grey = []
    edge_list_red = []
    edge_list_blue = []
    edge_list_purple = []
    for i in range(len(nodes)):
        n1 = nodes[i]
        for j in range(i+1,len(nodes)):
            n2 = nodes[j]
            if n1 == -n2:
                edge_list_grey.append((n1,n2))
            elif G2 is not None:
                if G.has_edge(n1,n2) and G2.has_edge(n1,n2):
                    edge_list_purple.append((n1,n2))
                elif G.has_edge(n1,n2):
                    edge_list_red.append((n1,n2))
                elif G2.has_edge(n1,n2):
                    edge_list_blue.append((n1,n2))                        
            elif G.has_edge(n1,n2):
                edge_list_red.append((n1,n2))
    nx.draw_networkx_edges(G, pos=pos,edgelist=edge_list_grey,edge_color='grey')
    nx.draw_networkx_edges(G, pos=pos,edgelist=edge_list_purple,edge_color='purple')
    if red_blue:
        nx.draw_networkx_edges(G, pos=pos,edgelist=edge_list_blue,edge_color='blue')
        nx.draw_networkx_edges(G, pos=pos,edgelist=edge_list_red,edge_color='red')
    else:
        nx.draw_networkx_edges(G, pos=pos,edgelist=edge_list_blue,edge_color='red')
        nx.draw_networkx_edges(G, pos=pos,edgelist=edge_list_red,edge_color='blue')
        
def show_combined(Gcombined,show_grey=True):
    red_edges = get_color_edges_combined(Gcombined,color="red")
    blue_edges = get_color_edges_combined(Gcombined,color="blue")
    purple_edges = get_color_edges_combined(Gcombined,color="purple")
    if show_grey:
        grey_edges = get_color_edges_combined(Gcombined,color="grey")
    pos = nx.circular_layout(Gcombined)
    nx.draw_networkx_nodes(Gcombined, pos=pos,node_size=600,node_color='w')
    nx.draw_networkx_labels(Gcombined, pos=pos)
            
    nx.draw_networkx_edges(Gcombined, pos=pos,edgelist=blue_edges,edge_color='blue')
    nx.draw_networkx_edges(Gcombined, pos=pos,edgelist=red_edges,edge_color='red') 
    nx.draw_networkx_edges(Gcombined, pos=pos,edgelist=purple_edges,edge_color='purple')

    if show_grey:
        nx.draw_networkx_edges(Gcombined, pos=pos,edgelist=grey_edges,edge_color='grey') 
    
def get_color_edges_combined(Gcombined,color="red"):
    color_edges = []
    df = pd.DataFrame(nx.adjacency_matrix(Gcombined).todense(),index=Gcombined.nodes(),columns=Gcombined.nodes())
    for i in range(len(df)):
        for j in range(len(df)):
            if df.iloc[i,j] == 1:
                data = Gcombined.get_edge_data(df.index[i],df.columns[j])
                if data['color'] == color:
                    color_edges.append((df.index[i],df.columns[j]))
    return color_edges

def genome_to_graph(genome):
    G = nx.Graph()
    ## BEGIN SOLUTION
    for P in genome:
        for i in range(len(P)):
            G.add_node(-1*P.iloc[i])
            G.add_node(P.iloc[i])
            G.add_edge(P.iloc[i],-1*P.iloc[i],color='grey')
        for i in range(len(P)-1):
            G.add_edge(P.iloc[i],-P.iloc[i+1],color='red')
        G.add_edge(P.iloc[i+1],-P.iloc[0],color='red')
    ## END SOLUTION
    return G

def combine(G,G2):
    Gcombined = nx.Graph()
    ## BEGIN SOLUTION
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        n1 = nodes[i]
        for j in range(i+1,len(nodes)):
            n2 = nodes[j]
            if n1 != -n2:
                if G2.has_edge(n1,n2) and G.has_edge(n1,n2):
                    Gcombined.add_edge(n1,n2,color="purple")
                elif G.has_edge(n1,n2):
                    Gcombined.add_edge(n1,n2,color="red")
                elif G2.has_edge(n1,n2):
                    Gcombined.add_edge(n1,n2,color="blue")
            else:
                if G.has_edge(n1,n2):
                    Gcombined.add_edge(n1,n2,color='grey')
                else:
                    Gcombined.add_edge(n2,n1,color='grey')
    ## END SOLUTION
    return Gcombined

def cycles(G,G2):
    nalt_cycles = 0
    ## BEGIN SOLUTION
    Gcombined = combine(G,G2)
    for edge in list(Gcombined.edges()):
        if edge[0] == -edge[1]:
            Gcombined.remove_edge(edge[0],edge[1])
    sub_graphs = [Gcombined.subgraph(c).copy() for c in nx.connected_components(Gcombined)] #nx.connected_component_subgraphs(Gcombined)
    nalt_cycles = len(list(sub_graphs))
    ## END SOLUTION
    return nalt_cycles

def blocks(G):
    ## BEGIN SOLUTION
    return round(len(list(G.nodes()))/2)
    ## END SOLUTION
    return 0

def two_break_distance(G,G2):
    ## BEGIN SOLUTION
    return blocks(G) - cycles(G,G2)
    ## END SOLUTION
    # YOUR SOLUTION HERE
    return 0

import matplotlib.pyplot as plt

def print_from_graph(G):
    sub_graphs = [G.subgraph(c).copy() for c in nx.connected_components(G)] #nx.connected_component_subgraphs(Gcombined)
    all_to_print = []
    for sub_graph in sub_graphs:   
        if len(list(sub_graph.nodes())) == 2:
            cycle = list(sub_graph.edges())
        else:
            cycle = list(nx.find_cycle(sub_graph))
        to_print = []
        for i in range(0,len(cycle),2):
            to_print.append(cycle[i][1])
        all_to_print.append(to_print)
    print("".join([str(c) for c in all_to_print]))
    return set([tuple(c) for c in all_to_print])
    
def get_color(sub_graph,edge):
    data = sub_graph.get_edge_data(edge[0],edge[1])
    return data['color']

def red_blue_cycle_check(sub_graph,cycle):
    checked_cycle = None
    colors = []
    ## BEGIN SOLUTION
    last_color = get_color(sub_graph,cycle[0])
    if last_color == 'grey':
        return None,None
    colors.append(last_color)
    for edge in cycle[1:]:
        color = get_color(sub_graph,edge)
        colors.append(color)
        if color == 'grey':
            return None,None
        if color == last_color:
            return None,None
        last_color = color
    checked_cycle = cycle
    ## END SOLUTION
    return checked_cycle,colors

def two_break_on_genome_graph(G,i1,i2,i3,i4,color='red'):
    G.remove_edge(i1,i2)
    G.remove_edge(i3,i4)
    if i1 != -i4:
        G.add_edge(i1,i4,color=color)
    if i2 != -i3:
        G.add_edge(i2,i3,color=color)
        
def shortest_rearrangement_scenario(P,Q):
    G_P = genome_to_graph(P)
    G_Q = genome_to_graph(Q)
    distance = two_break_distance(G_P,G_Q)
    Gcombined = combine(G_P,G_Q)
    fig = plt.figure(figsize=(20, 20));
    steps = [print_from_graph(G_P)]
    c=1
    plt.subplot(distance+1, 2, c); c+=1
    show_combined(Gcombined,show_grey=True)
    plt.subplot(distance+1, 2, c); c+=1
    show(G_P)#,P_G=Gcombined)
    first = True
    ## BEGIN SOLUTION
    while True:
        sub_graphs = connected_component_subgraphs(Gcombined)
        for sub_graph in sub_graphs:
            red_blue_cycle = None
            for cycle in nx.simple_cycles(sub_graph.to_directed()):
                edge_cycle = []
                a = cycle[0]
                for b in cycle[1:]:
                    edge_cycle.append([a,b])
                    a = b
                edge_cycle.append([b,cycle[0]])
                red_blue_cycle,colors = red_blue_cycle_check(sub_graph,edge_cycle)
                if red_blue_cycle is not None:
                    break
            if red_blue_cycle is None:
                return steps
            ix = colors.index('blue')
            #if first:
            #    ix = 4
            #    print(red_blue_cycle)
            #    first = False
            i1,i2 = red_blue_cycle[ix-1]
            i3,i4 = red_blue_cycle[ix+1]
        Gcombined.remove_edge(i1,i2)
        Gcombined.remove_edge(i3,i4)
        if i1 != -i4:
            Gcombined.add_edge(i1,i4,color='red')
        if i2 != -i3:
            Gcombined.add_edge(i2,i3,color='purple')
        two_break_on_genome_graph(G_P,i1,i2,i3,i4)
        plt.subplot(distance+1, 2, c); c+=1
        show_combined(Gcombined,show_grey=False)
        plt.subplot(distance+1, 2, c); c+=1
        show(G_P)#,P_G=Gcombined)
        steps.append(print_from_graph(G_P))
    ## END SOLUTION
    return steps