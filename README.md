# polymer-graph-simulation
Graph simulation of polymerization processes, allowing large flexibility. 

This package simulates the order of which monomers bond in a polymerisation process. At each step, the reactivity of every individual functional groups is found, and the combined reaction probability is calculated for each functional group pair. The reaction probability of a pair is simply the product of the functional group reactivities. The probabilities are normalised and used as weights to select a pair as the next reaction. 

Setting up a system has three steps. Below is an example for a two component copolymer system.

- Defining monomer types and the system composition. For a two component system with monomer functionalities 2 and 3, with 300 and 200 monomers used respectively:
  - functionalities = [2, 3]
  - loads = [300, 200]
  
- Defining the chemistry through kinetic reaction rates for each possible functional group pair:
  - For each component type - component type pair, a matrix holds the kinetic rates
  - In our example we have 3 of such pairs (A-A, B-B, A-B). If a monomer A can only react with B, the A-A and B-B kinetic rate tables will be 2x2 and 3x3 matrices with all elements equal to zero. If say one functional group on B _and_ A is twice as reactive as the others, the A-B kinetic rate table will be
  [[1, 1, 2]
   [2, 2, 4]]
 where the functional groups indexed as second and third are the more reactive groups.
 - Alternatively, two pre-defined sets are availalbe ('uniform_one_component' and 'uniform_two_component') where all rates are 1, and only A-B reaction are allowed for the latter
 
 - Selecting a First Shell Substitution Effect model
   - The original model of Flory ('flory') assumes each _monomer_ have the same reaction probability at all steps, effectively meaning that the reactivity of functional groups on reacted monomers increase
   - In the ideal polymerisation model ('id'/'ideal') the functional group reactivity does not change as the monomer it resides on reacts. Monomer reactivity is proportional to the number of unreacted groups (actually, the sum of their effective kinetic rate)
   - In the FSSE model, the functional group reactivity is a function of monomer degree (that holds the group). By default the FSSE model is defined by a list of arrays. Assuming that monomer A behaves ideally, but monomer B experience FSSE with the group reactivity halved for each degree increased, the FSSE list is
   FSSE = [[1, 1, 0], [4, 2, 1, 0]]
where the last item gives the reactivity of groups on a fully reacted monomer (0). Here a functional group on a free B monomers is 4 times as reactive as a group on a B monomer where that is the only reamining unreacted group. It is possible to include time dependence for the FSSE factors, using two files per monomer type: one that holds the FSSE values at given conversions in each line, the other containing the conversion values themselves (see example)

To see how a system is set up in practice, see the examples folder.
