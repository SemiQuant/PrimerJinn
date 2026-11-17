import random
import requests
import primer3
import concurrent.futures
import statistics
def generate_random_DNA_sequence(length):
    """Generate a random DNA sequence of a specified length,
    ensuring that no character is repeated more than twice in a row."""
    
    bases = ['A', 'C', 'G', 'T']
    seq = ''
    
    for i in range(length):
        # generate a random base
        base = random.choice(bases)
        
        # ensure that the same character does not appear more than twice in a row
        if len(seq) >= 2 and seq[-1] == seq[-2] == base:
            # generate a new base that is different from the previous two
            while base == seq[-1]:
                base = random.choice(bases)
        
        seq += base
    
    return seq

def q5_melting_temp(seq1, seq2="", salt=0.5):
    url = "https://tmapi.neb.com/tm/q5/%s/%s/%s?&fmt=long" % (salt, seq1, seq2)
    response = requests.get(url)
    json_string = response.json()
    tm1 = json_string["data"]["p1"][0]["tm"]
    return tm1

def process_sequence(seq):
    q5_tmp = q5_melting_temp(seq)
    me_tmp = primer3.bindings.calcTm(seq, dv_conc=2, mv_conc=70, dna_conc=3300)
    return q5_tmp, me_tmp

# q5_all = []
# me_all = []

q5 = []
me = []

with concurrent.futures.ThreadPoolExecutor(max_workers = 100) as executor:
    futures = []
    for i in range(15, 30):
        for j in range(1000):
            seq = generate_random_DNA_sequence(i)
            futures.append(executor.submit(process_sequence, seq))
    for future in concurrent.futures.as_completed(futures):
        q5_tmp, me_tmp = future.result()
        q5.append(q5_tmp)
        me.append(me_tmp)
        print(len(q5), end='\r')
        

# q5_all.append(q5)
# me_all.append(me)
# len(q5_all)
# q5 = [item for sublist in q5_all for item in sublist]
# me = [item for sublist in me_all for item in sublist]


list1 = q5
list2 = me

diff = [list1[i] - list2[i] for i in range(len(list1))]
# diff = [round(list1[i]) - round(list2[i]) for i in range(len(list1))]

mean_diff = statistics.mean(diff)
sd_diff = statistics.stdev(diff)

print("Mean difference:", mean_diff)
print("SD difference:", sd_diff)
