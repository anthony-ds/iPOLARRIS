#!/bin/bash

echo

realpath() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

ipoldir="$(realpath )"

sts=$1
ens=$2
raddir="$(realpath $3)"
tempdir="$(realpath $4)"
doppdir=$5
declare -a simdirs=(${@:6})

IFS=',' read -r -a starg <<< "$sts"
IFS=',' read -r -a enarg <<< "$ens"

ptype=mp4

mkdir -p $raddir $tempdir

declare -a mpopts=( "mp06" "mp08" "mp10" "mp16" "mp51" )
declare -a mpnames=( "wsm6" "thom" "morr" "wdm6" "p3" )

declare -a stt=()
declare -a edt=()

for ((ii=0;ii<${#starg[@]};ii++)); do
    st=$(echo ${starg[ii]} | sed 's/[^0-9]*//g')
    stdate=${st:0:8}
    sttime=${st:8:4}
    en=$(echo ${enarg[ii]} | sed 's/[^0-9]*//g')
    endate=${en:0:8}
    entime=${en:8:4}
    stt+=(${stdate}_${sttime})
    edt+=(${endate}_${entime})
done

configdir=$ipoldir/configtxt/${stt[0]}_${edt[${#edt[@]}-1]}
mkdir -p $configdir

echo Determining data agency and type...
echo
sleep 3

if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
	data='wrf'
    mp=$(ls $raddir/* | head -n 1 | xargs basename | cut -d '_' -f2)
    for ((ii=0;ii<${#mpopts[@]};ii++)); do
       [[ "${mpopts[ii]}" = "$mp" ]] && break
    done
    inputfile=input_${mpnames[ii]}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
    configfile=config_${mpnames[ii]}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
    echo $data
else
    if [ -z $simdirs ]; then
        data='obs'
        inputfile=input_${data}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
        configfile=config_${data}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
        echo $data
    else
        inputfile=input_obs_${stt[0]}_${edt[${#edt[@]}-1]}.txt
        configfile=config_obs_${stt[0]}_${edt[${#edt[@]}-1]}.txt
        echo "obs vs wrf"
    fi
fi

station=$(basename $raddir)
path2specs=$(realpath $raddir/../../../radar_specs/${station}_specs.txt)
latcen=$(cat $path2specs | grep "latitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
loncen=$(cat $path2specs | grep "longitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
minelev=$(cat $path2specs | grep "min_elev_angle =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
maxelev=$(cat $path2specs | grep "max_elev_angle =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
band=$(cat $path2specs | grep "frequency_band =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
agency=$(cat $path2specs | grep "agency =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs | tr '[:upper:]' '[:lower:]')

echo
echo Selecting radar files for analysis in range ${stt[0]} to ${edt[${#edt[@]}-1]}...
echo
sleep 3

tfile=input_${stt[0]}_${edt[${#edt[@]}-1]}.txt
rm -f $configdir/$tfile

for ((ii=0;ii<${#starg[@]};ii++)); do
 
    if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
        first=${stt[ii]:0:4}-${stt[ii]:4:2}-${stt[ii]:6:2}_${stt[ii]:9:2}
        last=${edt[ii]:0:4}-${edt[ii]:4:2}-${edt[ii]:6:2}_${edt[ii]:9:2}
        fmt=%Y-%m-%d_%H
        len=13
    else
        first="${stt[ii]:0:8}_${stt[ii]:9:2}"
        last="${edt[ii]:0:8}_${edt[ii]:9:2}"
        fmt=%Y%m%d_%H
        len=11
    fi
    dats=$(echo ${first:0:$len})
    while [[ $(echo $first | tr -d '-' | tr -d '_' | tr -d ':') -lt $(echo $last | tr -d '-' | tr -d '_' | tr -d ':') ]]; do
        first=$(date -j -v +1H -f $fmt $first +$fmt)
        dats=$dats"\|"${first:0:$len}
    done

    for filepath in $(ls $raddir/* | grep "$dats" | sort); do
        file=$(basename $filepath)
        if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
            filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
        else
            filedt=$(echo $file | cut -d '_' -f2)$(echo $file | cut -d '_' -f3)
        fi
        if [ "$filedt" -ge "$(echo ${stt[ii]} | tr -d '_')00" ] && [ "$filedt" -lt "$(echo ${edt[ii]} | tr -d '_')00" ]; then 
            echo $filepath >> $configdir/$tfile
            echo $(basename $filepath)
        fi
        if [ "$filedt" -ge "$(echo ${edt[ii]} | tr -d '_')00" ]; then
            break
        fi
    done

done

mv $configdir/$tfile $configdir/$inputfile

echo
echo Selecting temperature files for analysis in range ${stt[0]} to ${edt[${#edt[@]}-1]}...
echo
sleep 3 

tfile=tmp_${stt[0]}_${edt[${#edt[@]}-1]}.txt
rm -f $configdir/$tfile

if [[ "$(ls $tempdir/* | sort | head -n 1 | xargs basename)" == "wrfout"* ]]; then
    
    mp=$(ls $tempdir/* | sort | head -n 1 | xargs basename | cut -d '_' -f2)
    for ((ii=0;ii<${#mpopts[@]};ii++)); do
       [[ "${mpopts[ii]}" = "$mp" ]] && break
    done
    mpname=$(echo "${mpnames[ii]}")
    tempsrc=$mpname
    tempfile=temp_${tempsrc}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
    snd_on='False'
    wrft_on='True'

    for ((ii=0;ii<${#starg[@]};ii++)); do
        
        first=${stt[ii]:0:4}-${stt[ii]:4:2}-${stt[ii]:6:2}_${stt[ii]:9:2}
        last=${edt[ii]:0:4}-${edt[ii]:4:2}-${edt[ii]:6:2}_${edt[ii]:9:2}
        fmt=%Y-%m-%d_%H
        len=13
        dats=$(echo ${first:0:$len})
        while [[ $(echo $first | tr -d '-' | tr -d '_' | tr -d ':') -lt $(echo $last | tr -d '-' | tr -d '_' | tr -d ':') ]]; do
            first=$(date -j -v +1H -f $fmt $first +$fmt)
            dats=$dats"\|"${first:0:$len}
        done

        for filepath in $(ls $tempdir/* | grep "$dats" | sort); do
            file=$(basename $filepath)
            filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
            if [ "$filedt" -ge "$(echo ${stt[ii]} | tr -d '_')00" ] && [ "$filedt" -le "$(echo ${edt[ii]} | tr -d '_')00" ]; then 
                echo $filepath >> $configdir/$tfile
                echo $(basename $filepath)
            fi
            if [ "$filedt" -gt "$(echo ${edt[ii]} | tr -d '_')00" ]; then
                break
            fi
        done

    done

else

    tempsrc='uwyo-'$(basename $tempdir | tr '[:upper:]' '[:lower:]')
    tempfile=temp_${tempsrc}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
    snd_on='True'
    wrft_on='False'

    declare -a sndfiles=()

    for ((ii=0;ii<${#starg[@]};ii++)); do
        
        declare -a tmpsndfiles=()

        for filepath in $(ls $tempdir/* | sort); do

            file=$(basename $filepath)
            filedt=$(echo $file | cut -d '_' -f2 | tr -d '-')$(echo $file | cut -d '_' -f3 | cut -d '.' -f1 | tr -d ':')
            sttstr=$(echo ${stt[ii]} | tr -d '_')00
            edtstr=$(echo ${edt[ii]} | tr -d '_')00
            
            hrstdiff=$(( $(date -jf %Y%m%d%H%M%S "$filedt" +%s)-$(date -jf %Y%m%d%H%M%S "$sttstr" +%s) ))
            hrendiff=$(( $(date -jf %Y%m%d%H%M%S "$filedt" +%s)-$(date -jf %Y%m%d%H%M%S "$edtstr" +%s) ))

            if [ ${hrstdiff#-} -le 21600 ]; then
                tmpsndfiles+=($filepath)
            fi
            if [ ${hrendiff#-} -le 21600 ]; then
                tmpsndfiles+=($filepath)
            fi
        
        done

        if [ ! ${#tmpsndfiles[@]} -eq 2 ]; then
            echo Insufficient soundings for study period. Check and pull more!
            echo
            exit 1
        fi

        sndfiles+=($(echo ${tmpsndfiles[@]} | tr " " "\n" | sort -u | tr "\n" " "))

    done
    for file in ${sndfiles[@]}; do
        echo $file
        echo $file >> $configdir/$tfile
    done

fi

mv $configdir/$tfile $configdir/$tempfile

if [ ! -z $doppdir ]; then

    echo
    echo Selecting dual-Doppler files at closest times to range $stt to $edt...
    echo
    sleep 3

    dd_on='True'
    doppfile=dopp_${data}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
    
    tfile=dopp.txt
    rm -f $configdir/$tfile

    for ((ii=0;ii<${#starg[@]};ii++)); do

        if [[ "$(ls $doppdir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
            first=${stt[ii]:0:4}-${stt[ii]:4:2}-${stt[ii]:6:2}_${stt[ii]:9:2}
            last=${edt[ii]:0:4}-${edt[ii]:4:2}-${edt[ii]:6:2}_${edt[ii]:9:2}
            fmt=%Y-%m-%d_%H
            len=13
        else
            first="${stt[ii]:0:8}_${stt[ii]:9:2}"
            last="${edt[ii]:0:8}_${edt[ii]:9:2}"
            fmt=%Y%m%d_%H
            len=11
        fi
        dats=$(echo ${first:0:$len})
        while [[ $(echo $first | tr -d '-' | tr -d '_' | tr -d ':') -lt $(echo $last | tr -d '-' | tr -d '_' | tr -d ':') ]]; do
            first=$(date -j -v +1H -f $fmt $first +$fmt)
            dats=$dats"\|"${first:0:$len}
        done

        for filepath in $(ls $doppdir/* | grep "$dats" | sort); do
            file=$(basename $filepath)  
            if [[ "$(ls $doppdir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
                filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
            else
                filedt=$(echo $file | cut -d '_' -f3)$(echo $file | cut -d '_' -f4)
            fi
            if [ "$filedt" -ge "$(echo ${stt[ii]} | tr -d '_')00" ] && [ "$filedt" -lt "$(echo ${edt[ii]} | tr -d '_')00" ]; then 
                echo $filepath >> $configdir/$tfile
                echo $(basename $filepath)
            fi
            if [ "$filedt" -ge "$(echo ${edt[ii]} | tr -d '_')00" ]; then
                break
            fi
        done

    done

    mv $configdir/$tfile $configdir/$doppfile

else
    dd_on='False'
fi

if [ -z $simdirs ]; then
    if [[ "$data" == "obs" ]]; then
        fold="obs"
    else
        fold=$mpname
    fi

else

    declare -a inputfiles2=()
    declare -a configfiles2=()
    declare -a allmps=()

    for simdir in ${simdirs[@]}; do
        
        echo
        echo Selecting wrfout files for analysis in range $stt to $edt...

        fold='obsvwrf'
        mp2=$(ls $simdir/* | head -n 1 | xargs basename | cut -d '_' -f2)
        for ((ii=0;ii<${#mpopts[@]};ii++)); do
           [[ "${mpopts[ii]}" = "$mp2" ]] && break
        done
        mpname2=$(echo ${mpnames[ii]})
        echo MP=$(echo $mpname2 | tr '[:lower:]' '[:upper:]')
        echo

        inputfile2=input_${mpname2}_${stt[0]}_${edt[${#edt[@]}-1]}.txt
        configfile2=config_${mpname2}_${stt[0]}_${edt[${#edt[@]}-1]}.txt

        inputfiles2+=($configdir/$inputfile2)
        configfiles2+=($configdir/$configfile2)
        allmps+=($mpname2)

        tfile=input_${stt[0]}_${edt[${#edt[@]}-1]}.txt

        for ((ii=0;ii<${#starg[@]};ii++)); do
            
            first=${stt[ii]:0:4}-${stt[ii]:4:2}-${stt[ii]:6:2}_${stt[ii]:9:2}
            last=${edt[ii]:0:4}-${edt[ii]:4:2}-${edt[ii]:6:2}_${edt[ii]:9:2}
            fmt=%Y-%m-%d_%H
            len=13
            dats=$(echo ${first:0:$len})
            while [[ $(echo $first | tr -d '-' | tr -d '_' | tr -d ':') -lt $(echo $last | tr -d '-' | tr -d '_' | tr -d ':') ]]; do
                first=$(date -j -v +1H -f $fmt $first +$fmt)
                dats=$dats"\|"${first:0:$len}
            done

            for filepath in $(ls $simdir/* | grep "$dats" | sort); do
                file=$(basename $filepath)
                filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
                if [ "$filedt" -ge "$(echo ${stt[ii]} | tr -d '_')00" ] && [ "$filedt" -lt "$(echo ${edt[ii]} | tr -d '_')00" ]; then 
                    echo $filepath >> $configdir/$tfile
                    echo $(basename $filepath)
                fi
                if [ "$filedt" -ge "$(echo ${edt[ii]} | tr -d '_')00" ]; then
                    break
                fi
            done

        done

        mv $configdir/$tfile $configdir/$inputfile2

    done

fi

outfigdir=outputfig/${fold}_temp${tempsrc}_${station}_${stt[0]}_${edt[${#edt[@]}-1]}
outrrdir=$(cd $raddir/../../ && pwd)/radar_rainrates/$station
mkdir -p $outfigdir $outrrdir

template=$ipoldir/ipol_config.txt

if [ -z $simdir ]; then

    echo
    echo Creating config file for iPOLARRIS...
    echo
    sleep 3

    cp $template $configdir/$configfile

    sed -i '' "s/^type ==.*/type == $data == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
    if [[ "$data" == "obs" ]]; then
        sed -i '' "s/.*data ==.*/data == $agency == # Radar data source if type = 'obs' or model microphysics scheme if type = 'wrf'/g" $configdir/$configfile
    else
        sed -i '' "s/.*data ==.*/data == $mpname == # Radar data source if type = 'obs' or model microphysics scheme if type = 'wrf'/g" $configdir/$configfile
    fi
    sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile
    sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo ${stt[0]} | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo ${edt[${#edt[@]}-1]} | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
    if [[ "$wrft_on" == "True" ]]; then
        sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
    elif [[ "$snd_on" == "True" ]]; then
        sed -i '' "s%.*sfiles ==.*%sfiles == '$configdir/$tempfile' == # Path to list of sounding files to read in%g" $configdir/$configfile
    fi
    if [[ "$dd_on" == "True" ]]; then
        sed -i '' "s%.*dfiles ==.*%dfiles == '$configdir/$doppfile' == # Path to list of dual-Doppler files to read in%g" $configdir/$configfile
    fi
    if [[ "$data" == "obs" ]]; then
        sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
    else
        sed -i '' "s/.*exper ==.*/exper == $station-$(echo $mpname | tr '[:lower:]' '[:upper:]') == # Radar location/g" $configdir/$configfile
    fi
    sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/band ==.*/band == $band == # Frequency band/g" $configdir/$configfile
    sed -i '' "s/mincosthresh ==.*/mincosthresh == $minelev == # Min. elevation angle threshold for cone of silence/g" $configdir/$configfile
    sed -i '' "s/maxcosthresh ==.*/maxcosthresh == $maxelev == # Max. elevation angle threshold for cone of silence/g" $configdir/$configfile
    sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile
    sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile
    sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
    sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
    sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

    echo Running iPOLARRIS...
    sleep 3

    python run_ipolarris.py $configdir/$configfile

else

    echo
    echo Creating OBS and SIM config files for iPOLARRIS...
    echo
    sleep 3

    echo OBS
    cp $template $configdir/$configfile

    sed -i '' "s/^type ==.*/type == obs == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
    sed -i '' "s/.*data ==.*/data == $agency == # Radar data source if type = 'obs' or model microphysics scheme if type = 'wrf'/g" $configdir/$configfile
    sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile
    sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo ${stt[0]} | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo ${edt[${#edt[@]}-1]} | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
    if [[ "$wrft_on" == "True" ]]; then
        sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
    elif [[ "$snd_on" == "True" ]]; then
        sed -i '' "s%.*sfiles ==.*%sfiles == '$configdir/$tempfile' == # Path to list of sounding files to read in%g" $configdir/$configfile
    fi
    if [[ "$dd_on" == "True" ]]; then
        sed -i '' "s%.*dfiles ==.*%dfiles == '$configdir/$doppfile' == # Path to list of dual-Doppler files to read in%g" $configdir/$configfile
    fi
    sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
    sed -i '' "s/band ==.*/band == $band == # Frequency band/g" $configdir/$configfile
    sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
    sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile
    sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile
    sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
    sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
    sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

    for ((ii=0;ii<${#allmps[@]};ii++)); do
        
        echo MP=$(echo ${allmps[ii]} | tr '[:lower:]' '[:upper:]')
        cp $template ${configfiles2[ii]}

        sed -i '' "s/^type ==.*/type == wrf == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" ${configfiles2[ii]}
        sed -i '' "s/.*data ==.*/data == $(echo ${allmps[ii]}) == # Radar data source if type = 'obs' or model microphysics scheme if type = 'wrf'/g" ${configfiles2[ii]}
        sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" ${configfiles2[ii]}
        sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo ${stt[0]} | tr '_' '-')' == # Start time of analysis of interest/g" ${configfiles2[ii]}
        sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo ${edt[${#edt[@]}-1]} | tr '_' '-')' == # End time of analysis of interest/g" ${configfiles2[ii]}
        sed -i '' "s%.*rfiles ==.*%rfiles == '${inputfiles2[ii]}' == # Path to list of radar files to read in%g" ${configfiles2[ii]}
        sed -i '' "s%.*wfiles ==.*%wfiles == '${inputfiles2[ii]}' == # Path to list of WRF temperature files to read in%g" ${configfiles2[ii]}
        if [[ "$dd_on" == "True" ]]; then
            sed -i '' "s%.*dfiles ==.*%dfiles == '${inputfiles2[ii]}' == # Path to list of dual-Doppler files to read in%g" ${configfiles2[ii]}
        fi
        sed -i '' "s/.*exper ==.*/exper == $station-$(echo ${allmps[ii]} | tr '[:lower:]' '[:upper:]') == # Radar location/g" ${configfiles2[ii]}
        sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" ${configfiles2[ii]}
        sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" ${configfiles2[ii]}
        sed -i '' "s/band ==.*/band == $band == # Frequency band/g" ${configfiles2[ii]}
        sed -i '' "s/mincosthresh ==.*/mincosthresh == $minelev == # Min. elevation angle threshold for cone of silence/g" ${configfiles2[ii]}
        sed -i '' "s/maxcosthresh ==.*/maxcosthresh == $maxelev == # Max. elevation angle threshold for cone of silence/g" ${configfiles2[ii]}
        sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" ${configfiles2[ii]}
        sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" ${configfiles2[ii]}
        sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" ${configfiles2[ii]}
        sed -i '' "s/.*snd_on ==.*/snd_on == False == # Sounding temperature on/g" ${configfiles2[ii]}
        sed -i '' "s/.*wrft_on ==.*/wrft_on == True == # WRF temperature on/g" ${configfiles2[ii]}

    done

    echo
    echo Running iPOLARRIS...
    sleep 3

    python run_ipolarris.py $configdir/$configfile $(printf "%s " ${configfiles2[@]})

fi
