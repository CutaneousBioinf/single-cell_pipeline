import json
import logging
import os
from run_SoupX import parallel_SoupX
from run_Doublet import parallel_Doublet
from run_Merge import parallel_Merge
from run_QC import qc
from run_Norm import normalization
from run_PCA import pca
from run_Harmony import harmony
from run_Clustering import clustering


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename='log/info.log', filemode='w')

def main():
    os.chdir(os.getcwd())

    with open('config/config.json') as f:
        config = json.load(f)

    logging.info('------------------------------------------------------')
    logging.info('-------------------- Cooking Soup --------------------')
    logging.info('------------------------------------------------------')
    parallel_SoupX(config)

    logging.info('------------------------------------------------------')
    logging.info('------------------ Removing Doublet ------------------')
    logging.info('------------------------------------------------------')
    parallel_Doublet(config)

    logging.info('------------------------------------------------------')
    logging.info('------------------- Merging adatas -------------------')
    logging.info('------------------------------------------------------')
    parallel_Merge(config)

    logging.info('------------------------------------------------------')
    logging.info('------------------------ QCing -----------------------')
    logging.info('------------------------------------------------------')
    qc(config)

    logging.info('-------------------------------------------------------')
    logging.info('--------------------- Normalizing ---------------------')
    logging.info('-------------------------------------------------------')
    normalization(config)

    logging.info('-------------------------------------------------------')
    logging.info('------------------------ PCA -------------------------')
    logging.info('-------------------------------------------------------')
    pca(config)

    logging.info('-------------------------------------------------------')
    logging.info('--------------------- Harmony ------------------------')
    logging.info('-------------------------------------------------------')
    harmony(config)

    logging.info('-------------------------------------------------------')
    logging.info('-------------------- Clustering -----------------------')
    logging.info('-------------------------------------------------------')
    clustering(config)

    logging.info('-------------------------------------------------------')
    logging.info('----------------------- Done! -------------------------')
    logging.info('-------------------------------------------------------')


if __name__ == '__main__':
    main()