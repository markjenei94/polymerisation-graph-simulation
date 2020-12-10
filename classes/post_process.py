import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import glob


class PostProcess:

    def __init__(self, loads, reaction_input_type, **input_specs):
        '''
        :param reaction_input_type: str(array) or str(file)
        :param input_file_specs: reactions_folder_path=, reaction_file_root=, reactions=(if input_type == 'array')
        '''
        self.loads = loads
        if not (reaction_input_type == 'array' or reaction_input_type == 'file'):
            raise ValueError("'reaction_input_type' must be either 'array' or 'file'")
        self.reaction_input_type = reaction_input_type
        self.input_specs = input_specs

        self.reactions = self._load_reactions()
        self.G = nx.Graph()
        self.topology_data = {}

        self.number = []
        self.weight = []

        self.cycles_number = []
        self.cycles_mean_length = []

    def cycles(self, calculate_every=1, exclude_zeros=False):
        cycles_number = []
        cycles_mean_length = []

        for reactions in list(self.reactions.values()):
            print('.', end='')
            self._reset_graph()
            cycles_number.append([])
            cycles_mean_length.append([])
            i = 0
            for reaction in reactions:
                node1 = (reaction[0], reaction[2])
                node2 = (reaction[1], reaction[3])
                self.G.add_edge(node1, node2)
                if i % calculate_every != 0:
                    i += 1
                    continue
                i += 1
                cycles = nx.minimum_cycle_basis(self.G)
                number = len(cycles)
                mean_length = 0
                if number > 0:
                    mean_length = sum(len(x) for x in cycles) / float(number)
                cycles_number[-1].append(number)
                cycles_mean_length[-1].append(mean_length)
        print('')

        self.cycles_number = np.array(cycles_number)
        self.cycles_mean_length = np.array(cycles_mean_length)

        # mean_number_ave = np.mean(number_ave_data, axis=0)
        # mean_weight_ave = np.mean(weight_ave_data, axis=0)

        mean_cycles_number = self._mean(cycles_number)
        mean_cycles_mean_length = self._mean(cycles_mean_length, exclude_zeros=True)

        self.topology_data["cycles_number"] = mean_cycles_number
        self.topology_data["cycles_mean_length"] = mean_cycles_mean_length

    def average_sizes(self, reduced=False):
        number_ave_data = []
        weight_ave_data = []

        for reactions in list(self.reactions.values()):
            print('.', end='')
            self._reset_graph()
            number_ave_data.append([])
            weight_ave_data.append([])
            for reaction in reactions:
                node1 = (reaction[0], reaction[2])
                node2 = (reaction[1], reaction[3])
                self.G.add_edge(node1, node2)
                connected_components = [len(c) for c in sorted(nx.connected_components(self.G), key=len, reverse=True)]
                if reduced:
                    if len(connected_components) > 1:
                        connected_components = connected_components[1:]
                    else:
                        connected_components = [0.1]
                n_ave = np.average(connected_components)
                w_ave = np.average(connected_components, weights=connected_components)
                number_ave_data[-1].append(n_ave)
                weight_ave_data[-1].append(w_ave)
        print('')


        self.number = np.array(number_ave_data)
        self.weight = np.array(weight_ave_data)

        #mean_number_ave = np.mean(number_ave_data, axis=0)
        #mean_weight_ave = np.mean(weight_ave_data, axis=0)

        mean_number_ave = self._mean(number_ave_data)
        mean_weight_ave = self._mean(weight_ave_data)

        self.topology_data["number_average_size"] = mean_number_ave
        self.topology_data["weight_average_size"] = mean_weight_ave
        self.topology_data["dispersion_index"] = mean_weight_ave / mean_number_ave

    def _mean(self, array, exclude_zeros=False):
        lens = [len(i) for i in array]
        arr = np.ma.empty((np.max(lens), len(array)))
        arr.mask = True
        for idx, l in enumerate(array):
            arr[:len(l), idx] = l
        if exclude_zeros:
            arr = np.ma.masked_values(arr, 0.0)
        return arr.mean(axis=-1)

    def _load_reactions(self):
        loaded_reactions = {}
        if self.reaction_input_type == 'array':
            for i in range(len(self.input_specs['reactions'])):
                loaded_reactions[i] = self.input_specs['reactions'][i]
        else:
            dir_path = self.input_specs['reactions_folder_path']
            root = self.input_specs['reaction_file_root']
            files = glob.glob(str(dir_path) + str(root) + '*')
            i = 0
            for f in files:
                reactions = np.loadtxt(f)
                loaded_reactions[i] = reactions
                i += 1

        return loaded_reactions

    def _reset_graph(self):
        G = nx.Graph()
        for i in range(len(self.loads)):
            for j in range(self.loads[i]):
                G.add_node((i, j))
        self.G = G
