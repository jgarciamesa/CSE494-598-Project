import subprocess, re

class HLA:
    def __init__(self, A_1, A_2, B_1, B_2, C_1, C_2, DQA1_1, DQA1_2, DQB1_1, DQB1_2, DRB1_1, DRB2_2):
        self.A = (A_1, A_2)
        self.B = (B_1, B_2)
        self.C = (C_1, C_2)
        self.DQA1 = (DQA1_1, DQA1_2)
        self.DQB1 = (DQB1_1, DQB1_2)
        self.DRB1 = (DRB1_1, DRB2_2)


class Simulation:
    def __init__(self, simulation_id, n, known_HLA: HLA, typed_HLA: HLA):
        self.simulation_id = simulation_id
        self.n = n
        self.known_HLA = known_HLA
        self.typed_HLA = typed_HLA

    def num_allele_group_matches_by_gene(self):
        return self.num_regex_matches_by_gene('\w+\*(\d+):')
    
    def num_protein_sequence_matches_by_gene(self):
        return self.num_regex_matches_by_gene('\w+\*(\d+:\d+)')

    def num_regex_matches_by_gene(self, regex):
        num_regex_matches_by_gene = {}
        for gene in ['A', 'B', 'C', 'DQA1', 'DQB1', 'DRB1']:
            max_num_regex_matches_by_gene = 0
            
            known_allele_group_1 = re.search(regex, getattr(self.known_HLA, gene)[0]).group(1)
            known_allele_group_2 = re.search(regex, getattr(self.known_HLA, gene)[1]).group(1)

            typed_allele_group_1 = re.search(regex, getattr(self.typed_HLA, gene)[0]).group(1)
            typed_allele_group_2 = re.search(regex, getattr(self.typed_HLA, gene)[1]).group(1)
            
            num_uncrossed_regex_matches_by_gene = 0
            num_uncrossed_regex_matches_by_gene += 1 if known_allele_group_1 == typed_allele_group_1 else 0
            num_uncrossed_regex_matches_by_gene += 1 if known_allele_group_2 == typed_allele_group_2 else 0
            
            num_crossed_regex_matches_by_gene = 0
            num_crossed_regex_matches_by_gene += 1 if known_allele_group_1 == typed_allele_group_2 else 0
            num_crossed_regex_matches_by_gene += 1 if known_allele_group_2 == typed_allele_group_1 else 0

            num_regex_matches = max([num_uncrossed_regex_matches_by_gene, num_crossed_regex_matches_by_gene])
            num_regex_matches_by_gene[gene] = num_regex_matches
        return num_regex_matches_by_gene


find_simulations_output = subprocess.run('find . -mindepth 1 -maxdepth 1 -type d -name simulations', capture_output=True, text=True, shell=True).stdout.strip()

if not find_simulations_output:
    print("Please cd to main project directory! Exiting")
    exit(1)

simulations_ls_output = subprocess.run('ls -d simulations/simulation_*/', capture_output=True, text=True, shell=True).stdout.strip()
simulations = []
for simulation_directory_path in simulations_ls_output.splitlines():
    try:
        print(f'Processing {simulation_directory_path}')
        simulation_id = int(re.search(r'simulations/simulation_(\d+)', simulation_directory_path).group(1))

        find_output = subprocess.run(f'find simulations/simulation_{simulation_id}/results/fasta -name mixed_*.txt', capture_output=True, text=True, shell=True).stdout.strip()
        n = int(re.search(r'mixed_(\d+)', find_output).group(1))

        with open(f'simulations/simulation_{simulation_id}/results/fasta/mixed_{n}.txt') as known_profile:
            for i in range(4):
                known_profile.readline() # discard mtDNA info
            known_alleles = [re.search(r'\b\w+\*[\d:]+', known_profile_line).group(0) for known_profile_line in known_profile]
        known_hla_profile = HLA(*known_alleles)

        with open(f'simulations/simulation_{simulation_id}/results/kourami_{n}.result') as kourami_result:
            typed_alleles = [re.search(r'(^\w+\*[:\d]+)[; \w]', kourami_result_line).group(1) for kourami_result_line in kourami_result]
        typed_hla_profile = HLA(*typed_alleles)

        simulation = Simulation(simulation_id, n, known_hla_profile, typed_hla_profile)
        simulations.append(simulation)
    except:
        print(f'Processing of {simulation_directory_path} failed')

simulations.sort(key=(lambda simulation: simulation.n))
with open('simulations/hla_tally.txt', 'w') as hla_tally:
    for n in {simulation.n for simulation in simulations}:
        simulations_with_n_individuals = list(filter(lambda simulation: n == simulation.n, simulations))
        hla_tally.write(f'For {n} individuals:\n')
        num_protein_sequence_matches = sum([sum(simulation_with_n_individuals.num_protein_sequence_matches_by_gene().values()) for simulation_with_n_individuals in simulations_with_n_individuals])
        num_possible_matches = len(simulations_with_n_individuals) * 12
        hla_tally.write(f'There were {num_protein_sequence_matches} protein sequence matches (i.e. GENE*XX:YY resolution) out of {num_possible_matches} possible\n')

        num_allele_group_matches = sum([sum(simulation_with_n_individuals.num_allele_group_matches_by_gene().values()) for simulation_with_n_individuals in simulations_with_n_individuals])

        hla_tally.write(f'There were {num_allele_group_matches} allele group matches (i.e. GENE*XX resolution), incl. the above protein sequence matches, out of {num_possible_matches} possible\n\n')


