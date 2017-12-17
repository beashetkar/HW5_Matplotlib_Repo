
# PyMaceutical HW

Date: 12/16/17

The objective of this homework is to analyze the data to show how four treatments (Capomulin, Infubinol, Ketapril, and Placebo) compare. The following tasks are required:

1) Create a scatter plot that shows how the tumor volume changes over time for each treatment.
2) Create a scatter plot that shows how the number of metastatic (cancer spreading) sites changes over time for each treatment.
3) Create a scatter plot that shows the number of mice still alive through the course of treatment (Survival Rate).
4) Creating a bar graph that compares the total % tumor volume change for each drug across the full 45 days.      



```python
# Import Dependencies 
# Numpy for calculations and matplotlib for charting

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.stats import sem

```


```python
# CSV files to load
clinical_csv = os.path.join('raw_data', 'clinicaltrial_data.csv')
mouse_csv = os.path.join('raw_data', 'mouse_drug_data.csv')

# Read the data files into pandas
clinical_df = pd.read_csv(clinical_csv)
mouse_df = pd.read_csv(mouse_csv)

```


```python
clinical_df.count()
```


```python
# Drop duplicate rows in clinical_df except for the first ocurrence.

new_clinical_df=clinical_df.drop_duplicates(keep="first")
```


```python
new_clinical_df.count()
```


```python
# Merge the data into a single dataset.

total_data=pd.merge(new_clinical_df, mouse_df, how='inner',on='Mouse ID')

```


```python
total_data.columns
```


```python
total_data.head()
```


```python
# Create a compound index

total_data_df = total_data.set_index(["Mouse ID","Timepoint"])   # , "Tumor Volume" need 3 columns to make a unique index
```


## PART I .- Create a scatter plot that shows how the tumor volume changes over time for each treatment.



```python
# Grab DataFrame rows where column has certain values

drugs = ["Capomulin", "Infubinol", "Ketapril","Placebo" ]
drugs_df = total_data[total_data.Drug.isin(drugs)]
drugs_df.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
      <th>Drug</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>b128</td>
      <td>0</td>
      <td>45.000000</td>
      <td>0</td>
      <td>Capomulin</td>
    </tr>
    <tr>
      <th>1</th>
      <td>b128</td>
      <td>5</td>
      <td>45.651331</td>
      <td>0</td>
      <td>Capomulin</td>
    </tr>
    <tr>
      <th>2</th>
      <td>b128</td>
      <td>10</td>
      <td>43.270852</td>
      <td>0</td>
      <td>Capomulin</td>
    </tr>
    <tr>
      <th>3</th>
      <td>b128</td>
      <td>15</td>
      <td>43.784893</td>
      <td>0</td>
      <td>Capomulin</td>
    </tr>
    <tr>
      <th>4</th>
      <td>b128</td>
      <td>20</td>
      <td>42.731552</td>
      <td>0</td>
      <td>Capomulin</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Pivot Table for average 
tumor_change_mean = drugs_df.pivot_table(index=['Timepoint'], columns=['Drug'], values='Tumor Volume (mm3)', aggfunc=np.mean)

print (tumor_change.columns)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-8-a5ac9d6cf8d6> in <module>()
          2 tumor_change_mean = drugs_df.pivot_table(index=['Timepoint'], columns=['Drug'], values='Tumor Volume (mm3)', aggfunc=np.mean)
          3 
    ----> 4 print (tumor_change.columns)
    

    NameError: name 'tumor_change' is not defined



```python
# Pivot Table for standard error 
tumor_change_sem = drugs_df.pivot_table(index=['Timepoint'], columns=['Drug'], values='Tumor Volume (mm3)', aggfunc=sem)

print (tumor_change_sem.columns)
```

    Index(['Capomulin', 'Infubinol', 'Ketapril', 'Placebo'], dtype='object', name='Drug')
    


```python
# Generate the Plot 

#Plotting graphs (with error bars)
plt.errorbar(tumor_change_mean.index, tumor_change_mean.loc[:,'Capomulin'], yerr= tumor_change_sem.loc[:,'Capomulin'],
             fmt="o", markersize= 8, linestyle='dashed', linewidth=1, color="b")
plt.errorbar(tumor_change_mean.index, tumor_change_mean.loc[:,'Infubinol'], yerr= tumor_change_sem.loc[:,'Infubinol'],
             fmt="^", markersize= 8, linestyle='dashed', linewidth=1, color="r")
plt.errorbar(tumor_change_mean.index, tumor_change_mean.loc[:,'Ketapril'], yerr= tumor_change_sem.loc[:,'Ketapril'],
             fmt="s", markersize= 8, linestyle='dashed', linewidth=1, color="g")
plt.errorbar(tumor_change_mean.index, tumor_change_mean.loc[:,'Placebo'], yerr= tumor_change_sem.loc[:,'Placebo'],
             fmt="d", markersize= 8, linestyle='dashed', linewidth=1, color="k")

plt.title("Tumor Response over Treatment")
plt.xlabel("Time (days)")
plt.ylabel("Tumor volume (mm3)")
plt.grid(True)
plt.figsize=(20, 10)

plt.xlim(0,drugs_df['Timepoint'].max()+1)

plt.ylim(30,drugs_df['Tumor Volume (mm3)'].max()+5)

# Include a legend in best location
plt.legend(loc="best", fancybox=True)

# Save the figure
plt.savefig('charts/PyMaceutical_fig1.png')

# Show the figure
plt.show()

```


![png](PyMaceutical_files/PyMaceutical_14_0.png)


### Conclusions: "Tumor Response over Time"  ( PyMaceutical_fig1.png )

Notice the decrease in tumor volume change over the 45 days treatment for the mice that were treated with Capomulin. On other hand, the other two drugs were unsuccessful; been Ketrapril drug performed worst than a Placebo.



## PART II .- Create a scatter plot that shows how the number of metastatic (cancer spreading) sites changes over time for each treatment.



```python
# Pivot Table for average 
metastatic_mean = drugs_df.pivot_table(index=['Timepoint'], columns=['Drug'], values='Metastatic Sites', aggfunc=np.mean)


# Pivot Table for standard error 
metastatic_sem = drugs_df.pivot_table(index=['Timepoint'], columns=['Drug'], values='Metastatic Sites', aggfunc=sem)

```


```python
# Generate the Plot 

plt.title("Metastatic Spread During Treatment")
plt.xlabel("Treatment Duration (days)")
plt.ylabel("Metastatic Sites")
plt.grid(True)
plt.figsize=(20,10)

# Getting the limits
plt.xlim(0,drugs_df['Timepoint'].max()+1)
plt.ylim(0,drugs_df['Metastatic Sites'].max())

#Plotting graphs (with error bars)
plt.errorbar(metastatic_mean.index, metastatic_mean.loc[:,'Capomulin'], yerr= metastatic_sem.loc[:,'Capomulin'],
             fmt='o',  markersize= 8, linestyle='dashed', linewidth=1, color="b")
plt.errorbar(metastatic_mean.index, metastatic_mean.loc[:,'Infubinol'], yerr= metastatic_sem.loc[:,'Infubinol'],
             fmt="^",  markersize= 8, linestyle='dashed', linewidth=1, color="r")
plt.errorbar(metastatic_mean.index, metastatic_mean.loc[:,'Ketapril'], yerr= metastatic_sem.loc[:,'Ketapril'],
             fmt="s",  markersize= 8, linestyle='dashed', linewidth=1, color="g")
plt.errorbar(metastatic_mean.index, metastatic_mean.loc[:,'Placebo'], yerr= metastatic_sem.loc[:,'Placebo'],
             fmt="d",  markersize= 8, linestyle='dashed', linewidth=1, color="k")

# Include a legend in best location
plt.legend(loc="best", fancybox=True)

# Save the figure
plt.savefig('charts/PyMaceutical_fig2.png')

# Show the figure
plt.show()
```


![png](PyMaceutical_files/PyMaceutical_18_0.png)


### Conclusions: "Metastatic Spread During Treatment"  ( PyMaceutical_fig2.png )

In this chart, The drugs Ketapril and Placebo caused high risk of metastatic spread during the 45 days treatment. The best medication was Capomulin having lower risk, and the next one was Infubinol.



## PART III .- Create a scatter plot that shows the number of mice still alive through the course of treatment (Survival Rate).



```python
# Generate the Plot 

plt.style.use('ggplot')
plt.figsize=(20, 12)
plt.title("Survival During Treatment")
plt.xlabel("Treatment Duration (days)")
plt.ylabel("Survival Rate %")
plt.grid(True)

x =[0,5,10,15,20,25,30,35,40,45]

# Count how many times each mice appears in our group

drugs={}

for d in drugs_df.Drug.unique():
    drugs[d] = drugs_df.loc[drugs_df['Drug']==d,:]
    
    # Calculate the Percentage of mice alive
    survival_rate = drugs[d].groupby(['Drug','Timepoint'])['Tumor Volume (mm3)'].count()
    survival_rate = 100 *survival_rate /25
    
    # Plotting the Survival Rate for each Treatment
    #plt.scatter(x, survival_rate, label=d, alpha=0.6, linestyle='dashed', linewidth=1 )
    plt.plot(x, survival_rate, label=d, marker='o', markersize=8, alpha=0.6, linestyle='dashed', linewidth=1)
    
# Include a legend in best location
plt.legend(loc="best", fancybox=True)
         
# Save the figure
plt.savefig('charts/PyMaceutical_fig3.png')

# Show the figure
plt.show()
```


![png](PyMaceutical_files/PyMaceutical_21_0.png)


### Conclusions: "Survival During Treatment"  ( PyMaceutical_fig3.png )

The survival rate for mice treated with Capomulin was above 80% over 45 days of treatment, by contrast, the infunibol drug was not effective because the percentage of survival rate was below 50% after 30 days. The company should think to stop producing it.



## PART IV .- Creating a bar graph that compares the total % tumor volume change for each drug across the full 45 days.            



```python
# Use autolabel function to label the percentages of changes
def autolabel(rects, ax):
    
    # Get y-axis height to calculate label position from.
    (y_bottom, y_top) = ax.get_ylim()
    y_height = y_top - y_bottom

    for rect in rects:
        height = int(rect.get_height())
        
        # Fraction of axis height taken up by this rectangle
        p_height = (height / y_height)

        label_position = height - (y_height * 0.05)
        ax.text(rect.get_x() + rect.get_width()/2., height,'%d%%' % int(height), 
                ha='center', va='bottom', color='black')

```


```python
# Take the pivot table information of part I and calculate the percentage change for each drug
 
tumor_percentage_df = ((tumor_change_mean.iloc[-1] - tumor_change_mean.iloc[0])/tumor_change_mean.iloc[0])*100
drugs = ["Capomulin", "Infubinol", "Ketapril", "Placebo"]


fig, ax =plt.subplots()
plt.figsize=(20, 12)

# Set x values from dataframe 
x_values = np.arange(len(tumor_percentage_df))

# Plot the results, if passing color= green, if  failing drugs color = red

bar_width = 0.75
tumor_Pass = ax.bar(x_values[0], tumor_percentage_df[0], bar_width, linewidth= 1, linestyle='solid', color='green', alpha=0.5 )
tumor_Fail = ax.bar(x_values[1:], tumor_percentage_df[1:], bar_width, linewidth= 1, linestyle='solid' ,color='red', alpha=0.5 )

# Add labels to the x and y axes

ax.set_title("Tumor Volume Change over 45 Day Treatment")
ax.set_xlabel("Drugs")
ax.set_ylabel("% Tumor Volume Change")
ax.set_xticks(x_values + bar_width /2 )
ax.set_xticklabels(("Capomulin", "Infubinol", "Ketapril", "Placebo"))

ax.set_ylim(-25, tumor_percentage_df[1:].max()+10 )
tick_locations = [value for value in x_values]
plt.xticks(tick_locations, drugs)

# Add a semi-transparent horizontal line at y = 0
plt.hlines(0, -1, len(tumor_percentage_df))

ax.grid(True, linestyle='--')

# Call functions to implement the autolabel

autolabel(tumor_Pass, ax)
autolabel(tumor_Fail, ax)

# Include a legend in best location
plt.legend(loc="best", fancybox=True)

# Save the figure
plt.savefig('charts/PyMaceutical_fig4.png')

# Show the figure
plt.show()
```


![png](PyMaceutical_files/PyMaceutical_25_0.png)


### Conclusions: "Tumor Volume Change Over 45 Day Treatment"  ( PyMaceutical_fig4.png )

The  % of Tumor Change over the  45 day  treatment was above of 50% for Ketapril and Placebo, that means the tumor kept growing not matter if the mice was taking medication. The best drug over all was Capomulin; it performed really good because stop the increase of Tumor.    

