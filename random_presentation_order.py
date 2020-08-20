import secrets
#from numpy.random import shuffle

contestants = ["Filippo", "Jake", "Neha"]

print(f"Today's contestants are: {contestants}")

#shuffle(contestants)
secure_random = secrets.SystemRandom()
mixed_contest = secure_random.sample(contestants, len(contestants))
print(f"After randomizing them, we have:")
for count,person in enumerate(contestants):
    lucky_winner = mixed_contest[count]
    print(f"  Presenter #{count+1}: {lucky_winner}")
