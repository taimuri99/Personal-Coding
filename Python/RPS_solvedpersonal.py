# Incorporate the random library
import random

# Print Title
print("Let's Play Rock Paper Scissors!")

# Specify the three options
options = ["r", "p", "s"]

# Computer Selection
computer_choice = random.choice(options)

# User Selection
user_choice = input("Make your Choice: (r)ock, (p)aper, (s)cissors? ")

# Run Conditionals
if user_choice == "r" and computer_choice == "r":
    print("Tie")
elif user_choice == "r" and computer_choice == "p":
    print("Lose")
elif user_choice == "r" and computer_choice == "s":
    print("Win")
elif user_choice == "p" and computer_choice == "r":
    print("Win")
elif user_choice == "p" and computer_choice == "p":
    print("Tie")
elif user_choice == "p" and computer_choice == "s":
    print("Lose")
elif user_choice == "s" and computer_choice == "r":
    print("Lose")
elif user_choice == "s" and computer_choice == "p":
    print("Win")
elif user_choice == "s" and computer_choice == "s":
    print("Tie")
else:
    print("Wrong Command")
list = [1,2,3,4]
print(list[-1])
