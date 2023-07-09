import numpy as np

I = np.matrix([[1,0], [0,1]])
X = np.matrix([[0,1],[1,0]])
Z = np.matrix([[1, 0],[0, -1]])
Y = np.matrix([[0,-1j],[1j,0]])

def modadd(state, Z_or_X, index):
    if state[Z_or_X][index]:
        state[Z_or_X][index]=0
    else:
        state[Z_or_X][index]=1
    return state

def evolution(state, size):
    Z_part, X_part = state
    new_state = [[0 for _ in range(size)],[0 for _ in range(size)]]
    # X to Z evolution
    for i in range(size):
        if X_part[i]:
            # Z
            new_state = modadd(new_state, 0, i)
    # Z evolution
    for i in range(size):
        if Z_part[i]:
            if i == 0: # Boundary
                # Z
                new_state = modadd(new_state, 0, i)
                new_state = modadd(new_state, 0, i+1)
                # X
                new_state = modadd(new_state, 1, i)
            elif i < size-1:
                # Z
                new_state = modadd(new_state, 0, i-1)
                new_state = modadd(new_state, 0, i)
                new_state = modadd(new_state, 0, i+1)
                # X
                new_state = modadd(new_state, 1, i)
            else:
                # Z
                new_state = modadd(new_state, 0, i-1)
                new_state = modadd(new_state, 0, i)
                # X
                new_state = modadd(new_state, 1, i) 
    return new_state

def CM(A, B):
    return A*B-B*A


state_0=[[1, 0, 0, 0],[0, 0, 0, 0]] # (Z, X)

size = len(state_0[0])

state = evolution(state_0, size)

all_states = [state_0]

while state != state_0:
    all_states.append(state)
    state = evolution(state, size)

print(len(all_states))

paulis = [I, Z, X, Y]
p_set = []
for p_state in all_states:
    p_list = []
    Z_part, X_part = p_state
    for i in range(size):
        if X_part[i] and Z_part[i]:
            p_list.append(3)
        elif X_part[i]:
            p_list.append(2)
        elif Z_part[i]:
            p_list.append(1)
        else:
            p_list.append(0)
    if size > 2:
        add_p = np.kron(paulis[p_list[0]], paulis[p_list[1]])
    for i in range (2, size):
        add_p = np.kron(add_p, paulis[p_list[i]])
    p_set.append(add_p)

k = 4 # Qubit Count
d = 2**k # Qubit, 2

ortho = []
new_ortho = []

for p_mat in p_set:
    new = p_mat.copy()
    for mat in ortho:
        new = new-np.trace(mat.getH()*p_mat)*mat
    if new.any():
        ortho.append(new/np.sqrt(np.trace(new.getH()*new)/d))
        
        
# Generate Commutations

more = True

while more:
    commutations = []
    non_zero = []
    for mat1 in new_ortho: # THIS ONE IS TOO COMPUTATIONALLY DIFFICULT!!!!!
        for mat2 in ortho:
            if (mat1 - mat2).any():
                commutations.append(CM(mat1, mat2))

    new_ortho = []
    
    for mat_c in commutations:
        if mat_c.any():
            non_zero.append(mat_c)

    count = 0
    for mat_c in non_zero:
        new = mat_c.copy()
        for mat in ortho:
            new = new-(np.trace(mat.transpose().conjugate()*mat_c)/d)*mat
        if new.any():
            valid_mat = new/(np.sqrt(np.trace(new.transpose().conjugate()*new)/d))
            ortho.append(valid_mat)
            new_ortho.append(valid_mat)
            count += 1
    if count == 0:
        more = False
print(len(ortho))
