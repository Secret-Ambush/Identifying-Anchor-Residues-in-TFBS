def calculate_consensus(pwm):
    """
    Calculate the consensus sequence for a given PWM.
    
    Args:
    pwm (list of lists): The PWM as a list of lists, where each inner list represents
                         probabilities of residues at that position.
    
    Returns:
    str: The consensus sequence.
    """
    consensus = ''
    residues = ['A', 'C', 'G', 'T']
    
    for position in pwm:
        max_index = position.index(max(position))
        consensus += residues[max_index]
    
    return consensus

# Example PWM
pwm = [
    [0.542686, 0.175411, 0.161627, 0.120276],
    [0.336816, 0.140952, 0.160738, 0.361494],
    [0.338373, 0.147843, 0.141618, 0.372165],
    [0.383726, 0.126945, 0.133837, 0.355491],
    [0.373499, 0.136283, 0.149622, 0.340596],
    [0.369275, 0.133615, 0.129835, 0.367274],
    [0.379280, 0.169186, 0.136505, 0.315029],
    [0.572477, 0.204091, 0.121387, 0.102045]
]

consensus_sequence = calculate_consensus(pwm)
print(f"Consensus sequence: {consensus_sequence}")
