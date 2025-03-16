import numpy as np

# Creiamo un array vuoto con shape (4, 5, 10) e tipo object
mat = np.empty((4, 5, 10), dtype=object)

# Popoliamo l'array: con probabilità 50% un intero, altrimenti un float
for i in range(4):
    for j in range(5):
        for k in range(10):
            if np.random.rand() < 0.5:
                mat[i, j, k] = np.random.randint(0, 100)
            else:
                mat[i, j, k] = np.random.uniform(0, 100)

print("Matrice (4, 5, 10) con numeri misti:")
print(mat)

# Parametri iniziali
target_score = 30.0
# Coefficienti iniziali per 5 CV (il primo, ad esempio, ha un bias maggiore)
coeffs = np.array([1.2, 1.0, 1.0, 1.0, 1.0])

# Delta fissi (ipotetici) per ciascun CV: l'idea è che alcuni CV debbano aumentare, altri diminuire
delta = np.array([5.0, 2.0, 4.0, 1.0, 3.0])


def compute_score(coeffs, delta):
    """Calcola lo score come prodotto scalare tra i coefficienti e i delta."""
    return np.dot(coeffs, delta)


def reward(score):
    """Definisce il reward come il negativo della differenza assoluta rispetto al target."""
    return -abs(score - target_score)


# Parametri per l'ottimizzazione (simile a una semplice strategia di hill-climbing)
num_episodes = 100  # numero di iterazioni
step_size = 0.1  # variazione applicata ai coefficienti

# Valutazione iniziale
current_score = compute_score(coeffs, delta)
current_reward = reward(current_score)
print(f"Coefficients iniziali: {coeffs}")
print(f"Score iniziale: {current_score}, Reward: {current_reward}")

# Loop di "reinforcement learning" (approccio semplice)
for episode in range(num_episodes):
    # Copia dei coefficienti e scelta casuale di un indice da modificare
    new_coeffs = coeffs.copy()
    idx = np.random.randint(0, len(coeffs))
    # Modifica casuale: aumenta o diminuisce il coefficiente
    new_coeffs[idx] += np.random.choice([-step_size, step_size])

    # Calcola nuovo score e reward
    new_score = compute_score(new_coeffs, delta)
    new_reward = reward(new_score)

    # Se il reward migliora, accetta il nuovo set di coefficienti
    if new_reward > current_reward:
        coeffs = new_coeffs
        current_score = new_score
        current_reward = new_reward

print("\nCoefficients appresi:")
print(coeffs)
print(f"Score finale: {current_score}, Reward finale: {current_reward}")
