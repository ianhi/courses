ffmpeg -framerate 25 -i render_prep/force_shan_chen%06d.png -vf drawtext="fontfile=./OpenSans-Regular.ttf: text='Single Belt': fontcolor=white: fontsize=24: box=1: boxcolor=black@0.5: boxborderw=5: x=(w-text_w)/2: y=text_h/2" shan_chen_1_belt.mp4
ffmpeg -framerate 25 -i render_prep/force_shan_chen_double_belt%06d.png -vf drawtext="fontfile=./OpenSans-Regular.ttf: text='Double Belt': fontcolor=white: fontsize=24: box=1: boxcolor=black@0.5: boxborderw=5: x=(w-text_w)/2: y=text_h/2" shan_chen_2_belt.mp4
yes "y" | ffmpeg \
  -i shan_chen_1_belt.mp4 \
  -i shan_chen_2_belt.mp4 \
  -filter_complex '[0:v]pad=iw*2:ih[int];[int][1:v]overlay=W/2:0[vid]' \
  -map [vid] \
  -c:v libx264 \
  -crf 23 \
  compare_1_2_belt.mp4
