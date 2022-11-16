# Helper script containing a compressed screenign class used for deconvolution
import sys
import pandas as pd
import pickle as pl

class CompressedScreen:
	def __init__(self, metadata,design_matrix,features):

		"""
		metadata: pd.DataFrame, index is unique ids for samples columns are metadata
		design_matrix:  pd.DataFrame unique ids by perturbations
		features_matrix : pd.DataFrame, index is unique ids for samples columns are assay features
		"""
		self.metadata = metadata
		self.design_matrix = design_matrix
		self.features = features
		self.decon_result=None

	def subset(self, metadata_subset_columns, desired_values):
		"""
		Subset the metadata based on iterable metadata_subset_columns and 
		the corresponding desired values for each column.
		Only supports categorical columns
		------------------------------------------------------------
		Inputs:
			metadata_subset_columns: 	an iterable (e.g. list) of columns to subset on
			desired_values: 			an iterable of specific values for each column in 
										metadata_subset_columns to select.

		Return: 						a new object that is a subset of the original 
		"""
		metadata=self.metadata
		for i, meta_column in enumerate(metadata_subset_columns):
			metadata=metadata.loc[metadata[meta_column]==desired_values[i],:]
		
		design_matrix=self.design_matrix.loc[metadata.index, :]
		features=self.features.loc[metadata.index, :]
		subset_obj=CompressedScreen(metadata, design_matrix, features)
		return subset_obj
	
	def merge(self, other):
		"""
		Input:
			other: 			another CompressedScreen obj

		Return:
			merged_obj: 	a object that combined self and other
		"""

    	# make sure the unique ids don't collide
		unique_ids_1=self.metadata.index
		unique_ids_2=other.metadata.index
		assert len(set(unique_ids_1).intersection(set(unique_ids_2)))!=0, 'Unique IDs collide! Not allowed!'
		
		metadata=self.metadata.append(other.metadata)
		design_matrix=self.design_matrix.append(other.design_matrix)
		features=self.features.append(other.features)
		#append automatically writes NA to new columns
		#for the quantitative measures and design matrix, fill with zero
		features=features.fillna(0)
		design_matrix=design_matrix.fillna(0)
		merged_obj=CompressedScreen(metadata, design_matrix, features)
		return merged_obj

	def export(self, out_fn):
		"""
		Exports the object to pickle format
		"""
		pl.dump(self, open(out_fn, 'wb'))
	
	def save_deconvolution(self, results, type_):
		"""
		Saves the original deconvolution result to the class object

		Inputs:
			results: 		deconvolution result (e.g. matrix)
			type: 			a string describing what kinds of result (e.g. TSNE) or method  used
							to generate ('linear'), etc
		"""
		if self.decon_result==None:
			self.decon_result=dict()
		self.decon_result[type_]=results
	def get_deconvolution(self, type_):
		assert type_ in self.decon_result, 'Deconvolution type not found in object!'
		return self.decon_result[type_]
