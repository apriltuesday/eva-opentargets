# How to update workflow diagram

Workflow diagram is contained in [workflow.mmd](workflow.mmd). Its definition is written in [Mermaid](https://github.com/mermaid-js/mermaid).

## Regenerate PNG in command line

In order to regenerate the workflow diagram from source, use the commands:

```bash
# Install Mermaid CLI
npm install @mermaid-js/mermaid-cli

# Generate the diagram
./node_modules/.bin/mmdc -i workflow.mmd -s 3 -o workflow.png
```

## Live editor

If you want to adjust the diagram and see the result applied live, use [Mermaid live editor](https://mermaid-js.github.io/mermaid-live-editor/).
