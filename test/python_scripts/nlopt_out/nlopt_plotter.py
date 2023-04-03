import pandas as pd
import matplotlib.pyplot as plt


pathTabulardata = "/Users/tdinelli/Documents/GitHub/OptiSMOKE_toolbox/test/nlopt/test-batch/Output/optimization.out"
fun_type = "cm"

file = pd.read_csv(pathTabulardata, sep = '\s+', engine='python')


if(fun_type == "norm"):
	val = 1.000000e+07 # Furst
elif(fun_type == "cm"):	
	val = 1 # Narrow


obj = []
pen = []
ev = []
numpen = 0

for i in range(len(file['Obj-Function(19)'])):
	if file['Obj-Function(19)'][i] == val:
		numpen += 1
	if i == 0:
		obj.append(file['Obj-Function(19)'][i])
		pen.append(numpen)
		ev.append(i)
	elif file['Obj-Function(19)'][i] < obj[-1]:
		obj.append(file['Obj-Function(19)'][i])
		pen.append(numpen)
		ev.append(i)
	else:
		pass

fig, ax = plt.subplots()

ax.plot(ev, obj, 'g')
ax.set(xlabel = 'Evaluation number', ylabel = 'Obj function (1-CMscore)')
#ax.set_title(title)

# RECAP

print("=======================================")
#print(" ", title, "\n")
print("     Number of evaluation: ", len(file))
print("     First CM: ",1-obj[0])
print("     Last CM:  ",1-obj[-1])
print("=======================================")

plt.tight_layout()
plt.show()
