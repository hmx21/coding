import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class ReverseTuringMachine:
    def __init__(self, tape):
        # Initialize the Turing machine with the given tape
        self.tape = list(tape)
        self.head = 0  # Head starts at the beginning of the tape
        self.state = 'q0'  # Initial state
        self.states = {'q0', 'qf', 'r0', 'r1', 'l0', 'l1', 'l0\'', 'l1\''}  # Set of states
        self.alphabet = {' ', '0', '1', '*', '0\'', '1\''}  # Set of tape symbols
        # Transition function: (current_state, current_symbol) -> (new_state, new_symbol, direction)
        self.transitions = {
            # Beginning of the work
            ('q0', '0'): ('r0', '*', 1),  # Replace '0' with '*' and move right
            ('q0', '1'): ('r1', '*', 1),  # Replace '1' with '*' and move right
            ('q0', ' '): ('q0', ' ', -1),  # Replace '1' with '*' and move right

            # Transfer to the right
            ('r0', '0'): ('r0', '0', 1),  # Move right over '0'
            ('r0', '1'): ('r0', '1', 1),  # Move right over '1'
            ('r1', '0'): ('r1', '0', 1),  # Move right over '0'
            ('r1', '1'): ('r1', '1', 1),  # Move right over '1'

            # Change direction to the left
            ('r0', ' '): ('l0\'', ' ', -1),  # Move left if blank
            ('r1', ' '): ('l1\'', ' ', -1),  # Move left if blank
            ('r0', '0\''): ('l0\'', '0', -1),  # Move left over '0\''
            ('r1', '0\''): ('l1\'', '0', -1),  # Move left over '0\''
            ('r0', '1\''): ('l0\'', '1', -1),  # Move left over '1\''
            ('r1', '1\''): ('l1\'', '1', -1),  # Move left over '1\''

            # Place in the left adjacent cell
            ('l0\'', '0'): ('l0', '0\'', -1),  # Replace '0' with '0\'' and move left
            ('l1\'', '0'): ('l0', '1\'', -1),  # Replace '0' with '1\'' and move left
            ('l0\'', '1'): ('l1', '0\'', -1),  # Replace '1' with '0\'' and move left
            ('l1\'', '1'): ('l1', '1\'', -1),  # Replace '1' with '1\'' and move left

            # Transfer to the left
            ('l0', '0'): ('l0', '0', -1),  # Move left over '0'
            ('l1', '0'): ('l1', '0', -1),  # Move left over '0'
            ('l0', '1'): ('l0', '1', -1),  # Move left over '1'
            ('l1', '1'): ('l1', '1', -1),  # Move left over '1'

            # Change direction from left to right
            ('l0', '*'): ('q0', '0', 1),  # Replace '*' with '0' and move right
            ('l1', '*'): ('q0', '1', 1),  # Replace '*' with '1' and move right

            # Completion of work for even length
            ('q0', '0\''): ('qf', '0', -1),  # Replace '0\'' with '0' and move left to final state
            ('q0', '1\''): ('qf', '1', -1),  # Replace '1\'' with '1' and move left to final state

            # Completion of work for odd length
            ('l0\'', '*'): ('qf', '0', -1),  # Replace '*' with '0' and move left to final state
            ('l1\'', '*'): ('qf', '1', -1),  # Replace '*' with '1' and move left to final state
        }
        self.history = []  # To store the history of states, tape, and head positions

    def step(self):
        # Perform one step of the Turing machine
        if self.state == 'qf':
            return False  # Halt if in final state
        if self.head < 0 or self.head >= len(self.tape):
            self.tape.append(' ')  # Extend the tape with blank if head is out of bounds
        symbol = self.tape[self.head]
        if (self.state, symbol) in self.transitions:
            new_state, new_symbol, direction = self.transitions[(self.state, symbol)]
            self.tape[self.head] = new_symbol  # Write new symbol on the tape
            self.state = new_state  # Transition to new state
            self.head += direction  # Move head left or right
            # Save the current state, tape, and head position
            self.history.append((self.state, list(self.tape), self.head))
            return True
        return False

    def run(self):
        # Run the Turing machine until it halts
        while self.step():
            pass
        return ''.join(self.tape).strip()  # Return the final tape content

# Example usage !!! input the number of 0's and 1's you want to reverse
input_string = "010111001"
rtm = ReverseTuringMachine(input_string)
output_string = rtm.run()

# Animation setup
fig, ax = plt.subplots()
ax.set_title('Reverse Turing Machine')
tape_text = ax.text(0.5, 0.5, '', ha='center', va='center', fontsize=12, family='monospace', color='black')
state_text = ax.text(0.5, 0.9, '', ha='center', va='center', fontsize=12)
input_text = ax.text(0.1, 0.1, f'Input: {input_string}', ha='left', va='center', fontsize=12)
step_text = ax.text(0.9, 0.1, '', ha='right', va='center', fontsize=12)

# Create directory for saving images
output_dir = 'turing_machine_steps'
os.makedirs(output_dir, exist_ok=True)

def update(frame):
    # Update the animation for each frame
    state, tape, head = rtm.history[frame]
    tape_str = ''.join(tape)
    if 0 <= head < len(tape_str):
        tape_str = tape_str[:head] + '[' + tape_str[head] + ']' + tape_str[head+1:]
    tape_text.set_text(tape_str)
    state_text.set_text(f'State: {state}')
    step_text.set_text(f'Step: {frame}')
    if state == 'qf':
        ani.event_source.stop()  # Stop the animation if in final state
    # Save the current frame as an image
    plt.savefig(os.path.join(output_dir, f'step_{frame:03d}.png'))
    return tape_text, state_text, step_text

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(rtm.history), interval=500, blit=True)

# Save the animation as a gif using Pillow
ani.save('reverse_turing_machine.gif', writer='pillow')

plt.show()

print(f"Input: {input_string}")
print(f"Output: {output_string}")
