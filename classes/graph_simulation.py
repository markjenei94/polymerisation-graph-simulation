import numpy as np
import warnings
import functools
import pandas as pd
import networkx as nx


class KineticRates:

    def __init__(self, functionalities: list, preset=False):
        self.funcs = np.array(functionalities)
        self.dictionary = {}
        self.initialise_kinetic_rates(preset)

    def initialise_kinetic_rates(self, preset):
        implemented_presets = ["one_component_uniform", "two_component_uniform"]
        if not preset:
            return
        elif (preset == "one_component_uniform"):
            self.add_kinetic_rate_table(0, 0, np.ones([self.funcs[0], self.funcs[0]]))
        elif (preset == "two_component_uniform"):
            self.add_kinetic_rate_table(0, 0, np.zeros([self.funcs[0], self.funcs[0]]))
            self.add_kinetic_rate_table(1, 1, np.zeros([self.funcs[1], self.funcs[1]]))
            self.add_kinetic_rate_table(0, 1, np.ones([self.funcs[0], self.funcs[1]]))
        else:
            raise ValueError(
                f"Preset not defined. Set it to 'False' or choose from {implemented_presets}")

    def add_kinetic_rate_table(self, monomer_index1: int, monomer_index2: int, kinetic_rate_table):
        func1 = self.funcs[monomer_index1]
        func2 = self.funcs[monomer_index2]
        if (func1, func2) != np.array(kinetic_rate_table).shape:
            raise ValueError("Dimensions of kinetic rate table incorrect. != (funcs[index1], funcs[index2])")
        if (monomer_index1, monomer_index2) in self.dictionary.keys():
            warnings.warn("Kinetic rate table for these indices already loaded")
        self.dictionary[(monomer_index1, monomer_index2)] = np.array(kinetic_rate_table)

        if monomer_index1 != monomer_index2:
            self.dictionary[(monomer_index2, monomer_index1)] = np.array(kinetic_rate_table).T


class System:

    def __init__(self, loads: list, kinetic_rates: KineticRates):
        self.loads = np.array(loads)
        self.number_of_types = len(loads)
        self.kinetic_rates = kinetic_rates
        self.degree_dict = self.initialise_degrees()
        self.check_kinetic_rates()
        # self.print_kinetic_rates()

    def check_kinetic_rates(self):
        n = len(self.loads)
        missing_rates = []
        rates = self.kinetic_rates.dictionary
        funcs = self.kinetic_rates.funcs
        for i in range(n):
            for j in range(n):
                index = (i, j)
                if index not in rates:
                    self.kinetic_rates.add_kinetic_rate_table(i, j, np.zeros([funcs[i], funcs[j]]))
                    missing_rates.append(index)
        for missing_rate in missing_rates:
            warnings.warn("Missing kinetic rates for monomer type pairs " + str(missing_rate) + " initialised to zero")

    def print_kinetic_rates(self):
        np.set_printoptions(precision=3)
        print("| Kinetic Rate Table |")
        for key in self.kinetic_rates.dictionary.keys():
            print(f"\nMonomers {key[0]} and {key[1]}")
            print(self.kinetic_rates.dictionary[key])
            print("\n-----------------")

    def initialise_degrees(self):
        degree_dict = {}
        index = 0
        for monomer_func in self.kinetic_rates.funcs:
            degree_dict[index] = np.ones([self.loads[index], monomer_func])
            index += 1
        return degree_dict


class GraphSimulation:

    def __init__(self, system: System, model, cycles_allowed=True):
        self.system = system
        self.model = model
        self.cycles_allowed = cycles_allowed

        self.reactions = []

        if not cycles_allowed:
            self.G = self._initialise_graph()

        if model == 'data':
            self.fsse_values_file_paths = {}
            self.fsse_convs_file_paths = {}

        self.c = np.linspace(0, 1, int(np.dot(self.system.loads, self.system.kinetic_rates.funcs) / 2 + 1))
        self.step_ = 0

        self.log = {}

    def run(self, target_conversion=1., save_file_loc=''):
        self.fsse_factors = self._initialise_fsse_factors()
        if self.model == 'data':
            if (len(self.fsse_values_file_paths) == 0 or len(self.fsse_convs_file_paths) == 0):
                raise ValueError(
                    "Set the paths for the FSSE values and conversions files using self.set_fsse_file_paths(type, f_fsse, f_conv)")

        target = np.dot(self.system.loads, self.system.kinetic_rates.funcs) / 2
        if type(target_conversion) == float:
            target *= target_conversion
        elif type(target_conversion) == list:
            target_conversion = np.array(target_conversion)
            target = target * target_conversion
        else:
            raise TypeError(f"'target_conversion' must be of type {float} or {list}")

        below_target = True
        c = [0] * len(self.system.loads)
        step = 0
        p = 0
        pp = 0
        while below_target:
            reaction = self.step()
            if reaction == []:
                break

            c[reaction[0]] += 1
            c[reaction[1]] += 1
            step += 1
            self.step_ += 1

            p += 1
            if type(target_conversion) == float:
                if sum(c) / 2 >= target:
                    below_target = False
                if p > pp:
                    print('.', end='')
                    pp += target * 0.1
            else:
                if ((c >= target).any()):
                    below_target = False
                if p > pp:
                    print('.', end='')
                    pp += sum(target) * 0.1
        print('')

        if len(save_file_loc) != 0:
            reactions = np.array(self.reactions)
            np.savetxt(save_file_loc, reactions.astype(int), fmt='%i')

    def step(self):
        found = False
        while not found:
            functional_group_reactivities = self._functional_group_reactivities()
            pair_reactivites = self._functional_group_pair_reactivites(functional_group_reactivities)

            self.log["functional_group_reactivites"] = functional_group_reactivities
            self.log["pair_reactivities"] = pair_reactivites

            reaction_pair = self._sample_pair_reactivity_tensor(pair_reactivites)
            if len(reaction_pair) > 0:
                found = True
                if (reaction_pair[0] == reaction_pair[1] and reaction_pair[2] == reaction_pair[3]):
                    found = False
                if not self.cycles_allowed:
                    monomer1 = (reaction_pair[0], reaction_pair[2])
                    monomer2 = (reaction_pair[1], reaction_pair[3])
                    if self._cycles_pair_allowed(monomer1, monomer2) == False:
                        found = False

        self.reactions.append(reaction_pair)
        self._update_degrees(reaction_pair)
        return reaction_pair

    def set_fsse_file_paths(self, monomer_type, path_fsse, path_conv):
        self.fsse_values_file_paths[monomer_type] = path_fsse
        self.fsse_convs_file_paths[monomer_type] = path_conv

    def _initialise_fsse_factors(self):
        fsse_factors = {}
        i = 0
        for monomer_type in self.system.degree_dict:
            if (self.model == 'ideal' or self.model == 'id'):
                func = self.system.kinetic_rates.funcs[monomer_type]
                fsse_factors[monomer_type] = np.array(list(np.ones(func + 1)))
            elif type(self.model) == list:
                fsse_factors[monomer_type] = np.array(self.model[i])
                i += 1
            elif self.model == 'random':
                func = self.system.kinetic_rates.funcs[monomer_type]
                fsse_factors[monomer_type] = np.array(range(1, func + 2))
                fsse_factors[monomer_type][-1] = 0
            elif self.model == 'data':
                if monomer_type in self.fsse_values_file_paths:
                    fsse_factors[monomer_type] = self._interpolated_fsse_values(monomer_type)
                else:
                    func = self.system.kinetic_rates.funcs[monomer_type]
                    fsse_factors[monomer_type] = np.array(list(np.ones(func + 1)))

        return fsse_factors

    def _interpolated_fsse_values(self, monomer_type):
        fsse = []
        fsse_data = np.loadtxt(self.fsse_values_file_paths[monomer_type])
        fsse_data_c = np.loadtxt(self.fsse_convs_file_paths[monomer_type])
        for i in range(len(fsse_data)):
            fsse.append(np.interp(self.c, fsse_data_c, fsse_data[i]))
        fsse.append(np.zeros(len(fsse[-1])))
        return np.array(fsse).T

    def _update_degrees(self, reaction_pair: list):
        fg1 = reaction_pair[::2]
        fg2 = reaction_pair[1::2]
        self.system.degree_dict[fg1[0]][fg1[1]][fg1[2]] = 0
        self.system.degree_dict[fg2[0]][fg2[1]][fg2[2]] = 0

    def _functional_group_reactivities(self):
        functional_group_reactivities = {}
        for monomer_type in self.system.degree_dict:
            functional_group_reactivities[monomer_type] = self._fsse_corrected_degrees(monomer_type)
        return functional_group_reactivities

    def _fsse_corrected_degrees(self, monomer_type: int):
        if (self.model == 'id' or self.model == 'ideal'):
            return self.system.degree_dict[monomer_type]

        monomer_degrees = np.array(np.sum(1 - np.array(self.system.degree_dict[monomer_type]), axis=1)).astype(int)
        if self.fsse_factors[monomer_type].ndim == 1:
            fsse_factor = np.array([self.fsse_factors[monomer_type]] * len(self.system.degree_dict[monomer_type]))
        else:
            fsse_factor = np.array([self.fsse_factors[monomer_type][self.step_]] * len(self.system.degree_dict[monomer_type]))
        fsse_vector = fsse_factor.T[monomer_degrees].diagonal()
        fsse_vector = np.array([fsse_vector]).T

        fsse_corrected_degrees = self.system.degree_dict[monomer_type] * fsse_vector

        return fsse_corrected_degrees  # every monomer's every functional group corrected with monomer degree dep. fsse

    def _functional_group_pair_reactivites(self, functional_group_reactivities):
        pair_reactivites = {}
        n = len(self.system.degree_dict)
        # Using that types are assigned [0,n[ integers
        for type1 in range(n):
            for type2 in range(type1, n):
                kinetic_rates = self.system.kinetic_rates.dictionary[(type1, type2)]
                functional_group_reactivity_tensor = np.einsum('ij,jk->ijk', functional_group_reactivities[type1],
                                                               kinetic_rates)
                # monomer1-fgroup1, fgroup1-fgroup2 -> monomer1-fgroup1-fgroup2
                pair_reactivity = np.einsum('ijk,lk->iljk', functional_group_reactivity_tensor,
                                            functional_group_reactivities[type2])
                # monomer1-fgroup1-fgroup2, monomer2-fgroup2 -> monomer1-monomer2-fgroup1-fgroup2

                pair_reactivites[(type1, type2)] = pair_reactivity

        pair_reactivites = np.array(list(pair_reactivites.values()), dtype=object)
        return pair_reactivites

    def _sample_pair_reactivity_tensor(self, pair_reactivities):

        monomer_type_dict = {}
        index = 0
        reaction_index = 0
        # Using that types are assigned [0,n[ integers
        n = len(self.system.degree_dict)
        for type1 in range(n):
            for type2 in range(type1, n):
                s = np.sum(pair_reactivities[index])
                if s > 0:
                    monomer_type_dict[reaction_index] = list([type1, type2])
                    reaction_index += 1
                index += 1

        probabilities = []
        for reaction_type in pair_reactivities:
            s = np.sum(reaction_type)
            if s > 0:
                probabilities.append(reaction_type / np.sum(reaction_type))
        probabilities = np.array(probabilities)
        probabilities_flat = np.array(probabilities).flatten().astype(float)

        if sum(probabilities_flat) == 0:
            return []
        flat_index = np.random.choice(len(probabilities_flat), 1, p=probabilities_flat)
        tensor_indices = self._get_tensor_indices(probabilities, flat_index)
        reaction_indices = monomer_type_dict[int(tensor_indices[0])] + list(tensor_indices[1:].astype(int))

        return reaction_indices

    def _cartesian_product_broadcasted(self, *arrays):
        """
        http://stackoverflow.com/a/11146645/190597 (senderle)
        """
        broadcastable = np.ix_(*arrays)
        broadcasted = np.broadcast_arrays(*broadcastable)
        dtype = np.result_type(*arrays)
        rows, cols = functools.reduce(np.multiply, broadcasted[0].shape), len(broadcasted)
        out = np.empty(rows * cols, dtype=dtype)
        start, end = 0, rows
        for a in broadcasted:
            out[start:end] = a.reshape(-1)
            start, end = end, end + rows
        return out.reshape(cols, rows).T

    def _get_tensor_indices(self, pair_reactivities, flat_index):
        shape = pair_reactivities.shape
        coords = self._cartesian_product_broadcasted(*[np.arange(s, dtype='int') for s in shape])
        df = pd.DataFrame(coords)
        df['A'] = pair_reactivities.flatten()
        df = df.to_numpy()[:, :-1]
        return df[flat_index][0]

    def _cycles_pair_allowed(self, monomer1: tuple, monomer2: tuple):
        conversion = float(len(self.reactions)) / (np.dot(self.system.loads, self.system.kinetic_rates.funcs) / 2. )
        average_functionality = np.dot(self.system.loads, self.system.kinetic_rates.funcs) / np.sum(self.system.loads)
        gel_conversion_flory = 1. / (average_functionality - 1)
        gel_conversion_flory = 0.35
        if conversion < gel_conversion_flory:
            if nx.has_path(self.G, monomer1, monomer2):
                return False
            self.G.add_edge(monomer1, monomer2)
            return True
        elif self._cycles_pair_in_gel(monomer1, monomer2):
            self.G.add_edge(monomer1, monomer2)
            return True
        elif nx.has_path(self.G, monomer1, monomer2):
            return False
        self.G.add_edge(monomer1, monomer2)
        return True

    def _cycles_pair_in_gel(self, monomer1: tuple, monomer2: tuple):
        largest_component = max(nx.connected_components(self.G), key=len)
        if (monomer1 in largest_component or monomer2 in largest_component):
            return True
        return False

    def _initialise_graph(self):
        G = nx.Graph()
        for i in range(len(self.system.loads)):
            for j in range(self.system.loads[i]):
                G.add_node((i, j))
        return G