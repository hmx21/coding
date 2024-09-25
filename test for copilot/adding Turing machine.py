import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class BinaryAdditionTuringMachine:
    def __init__(self, tape):
        self.tape = list(tape)
        self.head = 0
        self.state = 'q0'
        self.states = {'q0', 'qf', 'dec', 'inc', 'borrow', 'carry', 'clear', 'r', 'l'}
        self.alphabet = {' ', '0', '1', '+'}
        self.transitions = {
            # Move to the right end of the tape
            ('q0', '0'): ('q0', '0', 1),
            ('q0', '1'): ('q0', '1', 1),
            ('q0', '+'): ('q0', '+', 1),
            ('q0', ' '): ('dec', ' ', -1),

            # Decrement the right number
            ('dec', '0'): ('borrow', '0', -1),
            ('dec', '1'): ('l', '0', -1),
            ('borrow', '0'): ('borrow', '0', -1),
            ('borrow', '1'): ('r', '0', 1),
            ('borrow', '+'): ('clear', ' ', 1),
            ('r', '0'): ('r', '1', 1),
            ('r', ' '): ('l', ' ', -1),

            # Move left to find the first bit to increment
            ('l', '0'): ('l', '0', -1),
            ('l', '1'): ('l', '1', -1),
            ('l', '+'): ('inc', '+', -1),

            # Handle increment
            ('inc', '0'): ('q0', '1', 1),
            ('inc', '1'): ('add', '0', -1),
            ('add', '0'): ('q0', '1', 1),
            ('add', '1'): ('add', '0', -1),
            ('add', ' '): ('carry', ' ', 1),

            # Clear the tape
            ('clear', '0'): ('clear', ' ', 1),
            ('clear', ' '): ('qf', ' ', 1),
        }
        self.history = []

    def step(self):
        if self.state == 'qf':
            return False
        if self.head < 0 or self.head >= len(self.tape):
            self.tape.append(' ')
        symbol = self.tape[self.head]
        if self.state == 'carry':
            self.tape.insert(0, '1')
            self.state = 'q0'
            self.history.append((self.state, list(self.tape), self.head))
            return True
        if (self.state, symbol) in self.transitions:
            new_state, new_symbol, direction = self.transitions[(self.state, symbol)]
            self.tape[self.head] = new_symbol
            self.state = new_state
            self.head += direction
            self.history.append((self.state, list(self.tape), self.head))
            return True
        return False

    def run(self):
        while self.step():
            pass
        return ''.join(self.tape).strip()

# Example usage
input_string = "101+11"
btm = BinaryAdditionTuringMachine(input_string)
output_string = btm.run()

# Animation setup
fig, ax = plt.subplots()
ax.set_title('Binary Addition Turing Machine')
tape_text = ax.text(0.5, 0.5, '', ha='center', va='center', fontsize=12, family='monospace', color='black')
state_text = ax.text(0.5, 0.9, '', ha='center', va='center', fontsize=12)
input_text = ax.text(0.1, 0.1, f'Input: {input_string}', ha='left', va='center', fontsize=12)
step_text = ax.text(0.9, 0.1, '', ha='right', va='center', fontsize=12)

def update(frame):
    state, tape, head = btm.history[frame]
    tape_str = ''.join(tape)
    if 0 <= head < len(tape_str):
        tape_str = tape_str[:head] + '[' + tape_str[head] + ']' + tape_str[head+1:]
    tape_text.set_text(tape_str)
    state_text.set_text(f'State: {state}')
    step_text.set_text(f'Step: {frame}')
    return tape_text, state_text, step_text

ani = animation.FuncAnimation(fig, update, frames=len(btm.history), interval=500, blit=True)

ani.save('binary_addition_turing_machine.gif', writer='pillow')

plt.show()

print(f"Input: {input_string}")
print(f"Output: {output_string}")
