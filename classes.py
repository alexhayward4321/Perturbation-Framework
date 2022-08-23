import numpy as np
import pandas as pd


class flux_object:
    """Just does some operations for me """

    def __init__(self, df):
        if self.consistent(df, 'particle'):
            self.particle = df['particle'][0]
        else:
            raise Exception("An error has occurred")
        self.cell = df['cell'][0]
        self.score = df['score'][0]

    def get_mid_bins(self, df):
        ...

    def consistent(df, label):
        if all(df[label] == df[label][0]):
            return True
        else:
            return False
