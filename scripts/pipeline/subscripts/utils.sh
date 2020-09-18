# Elaspsed time
convertsecs() {
  # Simple function for time elapsed
  h=$(bc <<<"${SECONDS}/3600")
  m=$(bc <<<"(${SECONDS}%3600)/60")
  s=$(bc <<<"${SECONDS}%60")
  printf "%02dH:%02dM:%05.2fS\n" "${h}" "${m}" "${s}"
}

echo_subsection() {
  echo -e "\n======== ${*} ========\n"
}

echo_section() {
  echo -e "++++++++ ${*} ++++++++\n"
}

