echo $1
OUTPUT=$1
FRAMERATE=24
echo $OUTPUT
identify -format "%h" ${OUTPUT}/initial.svg >> __H.txt 
identify -format "%w" ${OUTPUT}/initial.svg >> __W.txt 
expr 2 \* \( $(grep . __H.txt) / 2 \) >> __H1.txt 
expr 2 \* \( $(grep . __W.txt) / 2 \) >> __W1.txt 
echo "$(grep . __W1.txt)!x$(grep . __H1.txt)!" >> __resize.txt 
mogrify -format jpg -resize $(grep . __resize.txt) ${OUTPUT}/s*.svg
rm -f __H*.txt __W*.txt __resize.txt 
ffmpeg -r ${FRAMERATE} -f image2 -i ${OUTPUT}/snapshot%08d.jpg -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none ${OUTPUT}/out.mp4
##
rm ${OUTPUT}/*.jpg
#rm ${OUTPUT}/*.svg
#rm ${OUTPUT}/output*
	
