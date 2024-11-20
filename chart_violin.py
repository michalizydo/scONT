# Generation of violin plots from .arrow files from NanoComp.
# Usage : chart_violin.py

def main():

	import glob
	files=[]
	samplen=[]
	count=1
	
	for i in sorted(glob.glob("T7_MSA-1*.arrow")):
		samplen.append("MSA1 T7-"+str(count))
		files.append(i)
		count+=1
	count=1
	
	for i in sorted(glob.glob("MSA-1-*.arrow")):
		samplen.append("MSA1 RBP-"+str(count))
		files.append(i)
		count+=1
	count=1
	
	for i in sorted(glob.glob("Control-*.arrow")):
		samplen.append("Control RBP-"+str(count))
		files.append(i)
		count+=1
	count=1
	
	for i in sorted(glob.glob("MSA-2-*.arrow")):
		samplen.append("MSA2 RBP-"+str(count))
		files.append(i)
		count+=1
	count=1
	
	for i in sorted(glob.glob("P_Control-*.arrow")):
		samplen.append("Control T7-"+str(count))
		files.append(i)
		count+=1
	count=1
	
	for i in sorted(glob.glob("P_MSA-2-*.arrow")):
		samplen.append("MSA2 T7-"+str(count))
		files.append(i)
		count+=1
	
	import pandas as pd
	import pyarrow as pa
	names=[]
	alldata=[]
	for i in range(len(files)): 

		with pa.OSFile(files[i], 'rb') as source:
			loaded_array = pa.ipc.open_file(source).read_all()
		data=loaded_array.to_pandas()

		data=data['lengths'].values.tolist()
		# remove 80% outlier
		maxd=max(data)

		count5k=0
		for j in data:
			if j>3000:
				count5k+=1
			names.append(samplen[i])
		print(files[i]+" 3k:"+str(count5k))
		alldata+=data
	dataF=pd.DataFrame({'Experiment':names,'Read length':alldata})

	print(dataF)
	import seaborn
	import matplotlib.pyplot as plt
	plt.figure(figsize=(28,8),dpi=1200)
	seaborn.set(style = 'whitegrid') 

	seaborn.set(font_scale=1.6)
	ax = plt.gca()
	ax.set_ylim([-800, 150000])

if __name__=="__main__":
	main()