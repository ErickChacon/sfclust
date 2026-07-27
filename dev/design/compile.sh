#!/bin/bash
# Compile dev/design/design.md to a self-contained design.html with
# rendered Mermaid diagrams. Run from the repo root:
#   docker exec toolbox_rinla_latest bash -c "cd /home/rstudio/documents/repositories/sfclust && bash dev/design/compile.sh"
set -euo pipefail

cd "$(dirname "$0")"

quarto render design.md --to html -M embed-resources:true

MERMAID_DIR=$(find /usr/lib/rstudio-server/bin/quarto/share/formats/html/mermaid -maxdepth 1 -name mermaid.js -exec dirname {} \;)

# Quarto's HTML format ships mermaid.js/mermaid-init.js but only wires them up
# for {mermaid} cells (which need the quarto-ext/diagram extension). Plain
# ```mermaid fenced code blocks render to <pre class="mermaid"><code>...
# We render those client-side instead by switching to the "mermaid-js" class
# that mermaid-init.js scans for, and inlining the two scripts.
sed -i 's|<pre class="mermaid"><code>|<pre class="mermaid-js">|g; s|</code></pre>|</pre>|g' design.html

{
  echo '<script>'
  cat "$MERMAID_DIR/mermaid.js"
  echo '</script>'
  echo '<script>'
  cat "$MERMAID_DIR/mermaid-init.js"
  echo '</script>'
  echo '</body></html>'
} > /tmp/mermaid_tail.html

sed -i '/<\/body><\/html>/d' design.html
cat /tmp/mermaid_tail.html >> design.html
rm /tmp/mermaid_tail.html

echo "Compiled dev/design/design.html"
