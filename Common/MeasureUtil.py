import ecopy as ep
import numpy as np
from Common.FileUtil import TableData


class Diversity(object):
    def __init__(self, data: TableData, col_id: str, col_count: str):
        """"""
        self.data = data
        self.col_id = col_id
        self.col_count = col_count
        self.count = np.asarray([[d[col_count] for d in data]])
        self.count_sum = np.sum(self.count, axis=1)[0]

    def richness(self):
        """
        see https://en.wikipedia.org/wiki/Species_richness
        :return: 
        """
        return len(self.data)

    def shannons_entropy(self):
        """
        see https://en.wikipedia.org/wiki/Entropy_(information_theory)
        :return: 
        """
        return ep.diversity(self.count, 'shannon', num_equiv=False)[0]

    def gini_simpson_index(self):
        """
        see https://en.wikipedia.org/wiki/Diversity_index#Gini%E2%80%93Simpson_index
        :return: 
        """
        return ep.diversity(self.count, 'gini-simpson', num_equiv=False)[0]

    def simpson_index(self):
        """
        see https://en.wikipedia.org/wiki/Diversity_index#Simpson_index
        :return: 
        """
        return ep.diversity(self.count, 'simpson', num_equiv=False)[0]

    def berger_parker_dominance(self):
        """
        see https://en.wikipedia.org/wiki/Diversity_index#Berger%E2%80%93Parker_index
        :return: 
        """
        return ep.diversity(self.count, 'dominance', num_equiv=False)[0]

    def morisitas_index(self, data2: TableData, col_id2:str=None, col_count2:str=None) -> float:
        """
        see https://en.wikipedia.org/wiki/Morisita%27s_overlap_index
        :param data2: 
        :param col_id2: 
        :param col_count2: 
        :return: 
        """
        if col_id2 is None:
            col_id2 = self.col_id
        if col_count2 is None:
            col_count2 = self.col_count
        count_sum2 = np.sum([d[col_count2] for d in data2])

        map_id_counts = {}
        for d in self.data:
            map_id_counts[d[self.col_id]] = [d[self.col_count] / self.count_sum, 0]
        for d in data2:
            if d[col_id2] in map_id_counts:
                map_id_counts[d[col_id2]][1] = d[col_count2] / count_sum2
            else:
                map_id_counts[d[col_id2]] = [0, d[col_count2] / count_sum2]

        s = 0
        for xy in map_id_counts.values():
            s += xy[0] * xy[1]

        return 2 * s / (self.simpson_index() + Diversity(data2, col_id2, col_count2).simpson_index())

