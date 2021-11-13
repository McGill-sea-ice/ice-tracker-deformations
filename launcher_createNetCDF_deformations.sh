SCRIPT_CHARGE_DOMAINES="/home/smco810/REQA_plateforme/trunk/chargerDomainesUtilitaires.ksh"
LISTE_UTILITAIRES="SPI r.date"
. ${SCRIPT_CHARGE_DOMAINES} "${LISTE_UTILITAIRES}"


source .venv/bin/activate

start_date=2020040100
end_date=2020043000


aDate=${start_date}

while [[ $aDate -le ${end_date} ]]; do

yyyy=`echo ${aDate} | cut -c 1-4`
mm=`echo ${aDate} | cut -c 5-6`
dd=`echo ${aDate} | cut -c 7-8`


sed -i "s/start_year.*/start_year = ${yyyy}/" ./src/SeaIceDeformation/namelist.ini
sed -i "s/start_month.*/start_month = ${mm}/" ./src/SeaIceDeformation/namelist.ini
sed -i "s/start_day.*/start_day = ${dd}/" ./src/SeaIceDeformation/namelist.ini

aDate=$(r.date ${aDate} +24 | cut -c 1-10)

python src/SeaIceDeformation/main.py
echo "${aDate} is done"

done



