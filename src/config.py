'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Configuration script for data processing.
-----------------------------------------

Parameters:
overwrite -- (True, False) If overwrite is set to true, files that have already 
              been processed will be processed again and overwritten. 
              If overwrite is set to false, files that have aleady been processed 
              will not be processed during the program's execution.

'''

#-------------------- PARAMETERS -------------------------


overwrite = False


#----------------- INPUT .CSV FILES -----------------------

# List of raw .csv data paths to be processed by the main program
'''  
raw_csv_paths = ['/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200310211621_20200311184107_1.csv']

raw_csv_paths = [   "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200320010317_20200401010317_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200321163507_20200402163507_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200322153731_20200403153731_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200327041332_20200402041324_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200327055117_20200401055931_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200327122918_20200401123740_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200327181503_20200408181503_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200328013507_20200401010212_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200328013612_20200402014426_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200328013712_20200401010417_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200328031656_20200401033233_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200328040524_20200402041424_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200328180637_20200401182217_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301060745_20200302055937_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329021821_20200402014526_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329021921_20200402023422_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329121205_20200403122017_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329125928_20200401132356_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329210343_20200401105714_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329220112_20200410220112_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329220216_20200410220216_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200329220316_20200402213021_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330034653_20200401042046_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330034753_20200401033133_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330070806_20200402180128_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330084349_20200401082729_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330084449_20200401082829_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330111230_20200401105609_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330111534_20200402113945_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330125057_20200401123440_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330134116_20200401132456_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330134116_20200406133308_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330164856_20200402064141_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330191523_20200403193101_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330192620_20200401190957_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330192620_20200401195740_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330192720_20200401191057_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330192820_20200403194320_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330201200_20200401195540_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330201500_20200401195840_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330201600_20200401195940_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330210336_20200401204707_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330210336_20200401204812_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330210441_20200406205626_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330210541_20200411210541_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330214241_20200401105914_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330214341_20200401105814_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330214441_20200401105714_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331024742_20200401032833_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331024842_20200401024057_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331024942_20200401024202_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331025042_20200401024302_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331025142_20200401024402_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331042744_20200401041942_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331051822_20200401050931_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331060242_20200401032729_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331060242_20200402054622_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331060546_20200401050831_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331060646_20200401055931_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331065707_20200401064838_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331070007_20200402180128_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331074332_20200401064638_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331074432_20200401082529_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331074532_20200401073816_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331074632_20200401073916_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331074732_20200401074016_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331074832_20200403081306_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331074932_20200412074932_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331101635_20200405102532_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331101735_20200403210913_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331105838_20200401050727_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331105838_20200401064438_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331105943_20200401082429_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331110043_20200402104426_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331110725_20200401105914_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331115241_20200401123335_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331115346_20200402221000_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331115546_20200401114643_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331115646_20200401114743_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331124305_20200401123540_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331124405_20200401123640_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331141718_20200401131952_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331141822_20200401132056_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331141922_20200401132156_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331142022_20200402140400_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331142122_20200402140500_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331142222_20200403144651_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331150947_20200401150058_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331151051_20200401150158_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331151151_20200401150258_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331155713_20200401145958_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331160013_20200402163307_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331160113_20200401155257_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331164832_20200401164007_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331165036_20200401164107_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331172858_20200401064938_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331173335_20200401114143_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331173439_20200401145858_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331173639_20200401163907_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331173739_20200402172104_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331173939_20200404175626_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331174039_20200412174039_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331182820_20200401181917_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331182920_20200401182017_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331183020_20200409185419_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331190738_20200401082929_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331191348_20200401181612_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331191452_20200401181817_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331191552_20200401190853_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331195626_20200401194727_1.csv"]
    

raw_csv_paths = [   "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304013507_20200308010212_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304013611_20200309014426_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304013711_20200305021821_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304013711_20200308010416_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304031352_20200306025726_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304031456_20200306025830_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304031556_20200306025930_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304031656_20200306030030_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304040008_20200304053810_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304040319_20200305035606_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304040424_20200305035706_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304045237_20200305035501_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304045341_20200306043712_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304053810_20200305044059_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304053810_20200305053017_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304053915_20200304071753_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304054015_20200304103403_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304054115_20200305053346_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304063226_20200305062306_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304071753_20200304085555_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304071857_20200305062206_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304071957_20200304081006_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304072057_20200304170616_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304072157_20200305071435_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304072257_20200305183536_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304072357_20200305071635_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304081006_20200305071231_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304081111_20200304170616_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304081211_20200304175551_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304081311_20200304184335_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304085555_20200304103403_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304085659_20200304103508_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304085759_20200304103608_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304090319_20200305201252_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304103403_20200304121224_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304103403_20200305062006_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304103508_20200304121328_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304103608_20200305093921_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304104112_20200306111433_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304112803_20200304202422_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304112908_20200304121728_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304112908_20200305111943_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304113008_20200305112048_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304113108_20200305112248_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304121224_20200305044304_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304121328_20200305093821_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304121428_20200304170848_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304121728_20200305121005_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304121728_20200308123439_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304121928_20200304122028_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304130725_20200305125927_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304130925_20200305130027_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304135231_20200304153036_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304135435_20200305134728_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304135535_20200305134832_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304135635_20200305134932_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304153141_20200304170952_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304153241_20200305143542_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304153341_20200305152612_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304162348_20200304171252_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304162452_20200307160012_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304170616_20200305071335_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304170616_20200305080157_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304170848_20200304184714_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304170952_20200304184819_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304170952_20200305125523_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304171152_20200305175245_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304171252_20200304180232_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304171252_20200305161532_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304171352_20200306165729_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304175551_20200305080257_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180232_20200304185019_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180232_20200305175345_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180337_20200304185219_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180437_20200306174811_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180537_20200306183726_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180537_20200307174039_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180637_20200308182216_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304184131_20200305080457_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304184335_20200305080357_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304184714_20200305111805_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304184819_20200305143442_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304184919_20200304202718_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304185019_20200306183426_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304185219_20200306183526_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304185319_20200306183626_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304185419_20200306183726_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304194135_20200305193211_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304194240_20200305193311_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304194340_20200305202434_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304194340_20200306192720_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304194440_20200306192820_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304202422_20200305120901_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304202718_20200305175141_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304202822_20200305193111_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304202922_20200305202230_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304203022_20200304212002_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304210949_20200306205429_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304212002_20200306210336_1.csv",
                    "/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304212106_20200306210440_1.csv"]
'''

raw_csv_paths = ['/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200310211621_20200311184107_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331191348_20200401181612_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331060242_20200402054622_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330201200_20200401195540_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200330210336_20200401204707_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331182820_20200401181917_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331173639_20200401163907_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200331164832_20200401164007_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304171252_20200304180232_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304185019_20200306183426_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180232_20200305175345_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304180232_20200304185019_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304171152_20200305175245_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304202822_20200305193111_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304171352_20200306165729_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304203022_20200304212002_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304053810_20200305053017_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304053810_20200305044059_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304040008_20200304053810_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304194240_20200305193311_1.csv',
                  '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200304185319_20200306183626_1.csv']


def init():
    ''' () -> None

    Function that initializes the global lists that store .csv file paths for 
    every stage of data processing (processing, conversion and calculations), 
    using the global list of raw .csv data paths.

    Every raw .csv file (e.g. pairs_20200320010317_20200401010317_1.csv) is 
    associated to a processed (e.g. tri_20200320010317_20200401010317_1.csv), 
    converted (e.g. tri_gridCS_20200320010317_20200401010317_1.csv) and 
    calculated (e.g. calc_20200320010317_20200401010317_1.csv) .csv file.

    The .csv files are stored in the following tree structure, where 
    some_parent_folder can have any name :

    ├── 2021_SeaIceDeformations
    |    ...
    |    ├── data
    |    |   ├── 00_grid /etc
    |    |   ├── 01_raw
    |    |   |   └── some_parent_folder
    |    |   |       └── pairs_20200320010317_20200401010317_1.csv
    |    |   ├── 02_processed
    |    |   |   └── some_parent_folder
    |    |   |       └── tri_20200320010317_20200401010317_1.csv
    |    |   ├── 03_converted
    |    |   |   └── some_parent_folder
    |    |   |       └── tri_gridCS_20200320010317_20200401010317_1.csv
    |    |   └── 04_calculations
    |    |       └── some_parent_folder
    |    |           └── calc_20200320010317_20200401010317_1.csv
    ...
    
    '''

    from src.d00_utils.get_data_paths import get_processed_csv_path, get_converted_csv_path, get_calculations_csv_path

    # Initialize global variables as empty lists
    global processed_csv_paths, converted_csv_paths, calculations_csv_paths
    
    processed_csv_paths     = []
    converted_csv_paths     = []
    calculations_csv_paths  = []

    # Iterate through all raw .csv file paths
    for raw_csv_path in raw_csv_paths:
        # For each raw .csv file path, find the appropriate path for each subsequent stages of processing
        processed_csv_path, _    = get_processed_csv_path(raw_csv_path)     # path for processed .csv data file
        converted_csv_path, _    = get_converted_csv_path(raw_csv_path)     # path for converted .csv data file
        calculations_csv_path, _ = get_calculations_csv_path(raw_csv_path)  # path for calculated .csv data file

        # Add the .csv file paths to the file path lists
        processed_csv_paths.append(processed_csv_path)
        converted_csv_paths.append(converted_csv_path)
        calculations_csv_paths.append(calculations_csv_path)
