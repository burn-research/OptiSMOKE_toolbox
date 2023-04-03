import pandas as pd
import matplotlib.pyplot as plt
import json

config = json.load(open("tabdata_plotter.config"))
paths = config["PATH"]
names = config["NAMES"]
obj_fun_type = config["OBJ_FUN_TYPE"]

dimension_subplot = len(paths)

if dimension_subplot > 1:
	fig, axs = plt.subplots(dimension_subplot)
else:
	fig, axs = plt.subplots()

for k in range(len(paths)):

	pathTabulardata = paths[k]
	title = names[k]
	fun_type = obj_fun_type[k]

	file = pd.read_csv(pathTabulardata, sep = '\s+', engine='python')

	# print(file['obj_fn'])

	if(fun_type == "norm"):
		val = 1.000000e+07 # Furst
	elif(fun_type == "cm"):	
		val = 1 # Narrow

	# count = file['obj_fn'].value_counts()[val] # number_pen nel matlab

	obj = []
	pen = []
	ev = []
	numpen = 0

	for i in range(len(file['obj_fn'])):
		if file['obj_fn'][i] == val:
			numpen += 1
		if i == 0:
			obj.append(file['obj_fn'][i])
			pen.append(numpen)
			ev.append(i)
		elif file['obj_fn'][i] < obj[-1]:
			obj.append(file['obj_fn'][i])
			pen.append(numpen)
			ev.append(i)
		else:
			pass

	# print(ev)
	# print(obj)
	# Check se tutto funziona non dovrebbe causare problemi!
	# if numpen != count:
	# 	raise ValueError('There are errors in the above for cycle: \n The number of the calculated penalties must be equal to the number of reported penalties')

	if dimension_subplot > 1:
		axs[k].plot(ev, obj, 'g')
		axs[k].set(xlabel = 'Evaluation number', ylabel = 'Obj function (1-CMscore)')
		axs[k].set_title(title)
	else:
		axs.plot(ev, obj, 'g')
		axs.set(xlabel = 'Evaluation number', ylabel = 'Obj function (1-CMscore)')
		axs.set_title(title)

	# RECAP

	print("=======================================")
	print(" ", title, "\n")
	print("     Number of evaluation: ", len(file))
	print("     First CM: ",1-obj[0])
	print("     Last CM:  ",1-obj[-1])
	print("=======================================")

plt.tight_layout()
plt.show()
