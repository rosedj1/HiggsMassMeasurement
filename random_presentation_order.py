from numpy.random import shuffle

contestants = ["Filippo", "Jake", "Neha"]

shuffle(contestants)

print(f"Today's contestants are: {contestants}")
for count,person in enumerate(contestants):
    lucky_winner = contestants[count]
    print(f"  Presenter #{count+1}: {lucky_winner}")