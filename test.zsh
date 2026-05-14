#!/usr/bin/env zsh

# Find every Article_Reference.md and prepend a top-level header based on parent folder name.
# Example:
#   ./my_project_folder/Article_Reference.md
# becomes header:
#   # My Project Folder

set -euo pipefail

find . -type f -name 'Article_Reference.md' | while IFS= read -r file; do
  dir_name="${file:h:t}"

  # Replace underscores with spaces
  title="${dir_name//_/ }"

  # Convert to title case
  title=$(print -r -- "$title" | awk '{
    for (i = 1; i <= NF; i++) {
      $i = toupper(substr($i,1,1)) tolower(substr($i,2))
    }
    print
  }')

  tmp_file="$(mktemp)"

  {
    print -r -- "# $title"
    cat "$file"
  } > "$tmp_file"

  mv "$tmp_file" "$file"
  print "Updated: $file"
done