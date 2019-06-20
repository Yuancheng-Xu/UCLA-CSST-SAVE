#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')

import random
from pprint import pprint

from decision_tree_functions import decision_tree_algorithm, decision_tree_predictions
from helper_functions import train_test_split, calculate_accuracy


# # Load and Prepare Data

# #### Format of the data
# - last column of the data frame must contain the label and it must also be called "label"
# - there should be no missing values in the data frame

# In[7]:


df = pd.read_csv("winequality-red.csv")
df["label"] = df.quality
df = df.drop("quality", axis=1)

column_names = []
for column in df.columns:
    name = column.replace(" ", "_")
    column_names.append(name)
df.columns = column_names

df.head()


# In[14]:


# counts all the label and sort them
df.label.value_counts(normalize=True)


# In[6]:


wine_quality = df.label.value_counts(normalize=True)
wine_quality = wine_quality.sort_index()
wine_quality.plot(kind="bar")


# In[15]:


def transform_label(value):
    if value <= 5:
        return "bad"
    else:
        return "good"

df["label"] = df.label.apply(transform_label)


# In[8]:


wine_quality = df.label.value_counts(normalize=True)
wine_quality[["bad", "good"]].plot(kind="bar")
wine_quality


# In[19]:


random.seed(0)
train_df, test_df = train_test_split(df, test_size=0.2)


# # Random Forest

# In[22]:


# randomly select n_bootstrap data points with replacement (allows duplicates) 
def bootstrapping(train_df, n_bootstrap):
    bootstrap_indices = np.random.randint(low=0, high=len(train_df), size=n_bootstrap)
    df_bootstrapped = train_df.iloc[bootstrap_indices]
    
    return df_bootstrapped

# see decision_tree_functions.py for randomly select n_features features
def random_forest_algorithm(train_df, n_trees, n_bootstrap, n_features, dt_max_depth):
    forest = []
    for i in range(n_trees):
        df_bootstrapped = bootstrapping(train_df, n_bootstrap)
        tree = decision_tree_algorithm(df_bootstrapped, max_depth=dt_max_depth, random_subspace=n_features)
        forest.append(tree)
    
    return forest

def random_forest_predictions(test_df, forest):
    df_predictions = {}
    for i in range(len(forest)):
        column_name = "tree_{}".format(i)
        predictions = decision_tree_predictions(test_df, tree=forest[i])
        # keys and values
        df_predictions[column_name] = predictions
    
    # transform the dictionary into dataframe. 
    # rows:index of datapts; columns: prediction of each trees
    df_predictions = pd.DataFrame(df_predictions)
    # vote for most. Type of predicion is pandas.core.series.Series
    random_forest_predictions = df_predictions.mode(axis=1)[0]
    
    return random_forest_predictions


# In[24]:


forest = random_forest_algorithm(train_df, n_trees=20, n_bootstrap=400, n_features=5, dt_max_depth=4)
predictions = random_forest_predictions(test_df, forest)
accuracy = calculate_accuracy(predictions, test_df.label)

print("Accuracy = {}".format(accuracy))

