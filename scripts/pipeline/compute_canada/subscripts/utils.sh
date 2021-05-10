# Print text for subsection with date and runtime
echo_subsection() {
  h=$(bc <<<"${SECONDS}/3600")
  m=$(bc <<<"(${SECONDS}%3600)/60")
  s=$(bc <<<"${SECONDS}%60")
  echo -e "\n======== [$(date +%H:%M:%S)] ${*} [Runtime: ${h}h ${m}m ${s}s] ========\n"
}

# Print text for section with date and runtime
echo_section() {
  h=$(bc <<<"${SECONDS}/3600")
  m=$(bc <<<"(${SECONDS}%3600)/60")
  s=$(bc <<<"${SECONDS}%60")
  echo -e "\n++++++++ [$(date +%H:%M:%S)] ${*} [Runtime: ${h}h ${m}m ${s}s] ++++++++\n"
}
