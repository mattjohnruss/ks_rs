ffmpeg \
    -pattern_type glob \
    -i 'spatial_t=*_colour=m_i_factor.png' \
    -vf "drawtext=fontfile=Arial.ttf: text='%{frame_num}': start_number=0: x=2*lh: y=2*lh: fontcolor=black: fontsize=40: box=1: boxcolor=white: boxborderw=5" \
    -r 10 \
    -y m_i_factor.mp4

ffmpeg \
    -pattern_type glob \
    -i 'spatial_t=*_colour=j_phi_i_i_factor.png' \
    -vf "drawtext=fontfile=Arial.ttf: text='%{frame_num}': start_number=0: x=2*lh: y=2*lh: fontcolor=black: fontsize=40: box=1: boxcolor=white: boxborderw=5" \
    -r 10 \
    -y j_phi_i_i_factor.mp4

ffmpeg \
    -pattern_type glob \
    -i 'spatial_t=*_colour=t_j_phi_i_lag.png' \
    -vf "drawtext=fontfile=Arial.ttf: text='%{frame_num}': start_number=0: x=2*lh: y=2*lh: fontcolor=black: fontsize=40: box=1: boxcolor=white: boxborderw=5" \
    -r 10 \
    -y t_j_phi_i_lag.mp4
