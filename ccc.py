num = range(11)
lab = [chr(int("258"+str(i), 16)) if i > 0 else " " for i in num]
print(lab)