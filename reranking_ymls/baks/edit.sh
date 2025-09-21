for f in *; do
  [ -f "$f" ] && sed -i "/vdg_lib_dir: '\/wynton\/group\/degradolab\/skt\/docking\/databases\/frag_vdg_lib\/'/d" "$f"
done

