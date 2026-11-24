// Follow MkDocs Material's palette when embedded and the system palette when opened directly.
const schemePreference = matchMedia("(prefers-color-scheme: dark)");
const frame = window.frameElement;
const parentBody = frame?.ownerDocument.body;

function updateTheme() {
  const scheme = parentBody?.getAttribute("data-md-color-scheme");
  const dark = scheme ? scheme === "slate" : schemePreference.matches;
  document.documentElement.dataset.theme = dark ? "dark" : "light";
}

updateTheme();
schemePreference.addEventListener("change", updateTheme);
if (frame) {
  new MutationObserver(updateTheme).observe(parentBody, {
    attributes: true,
    attributeFilter: ["data-md-color-scheme"],
  });
  // Let readers scroll the documentation page, including when they expand Model details.
  new ResizeObserver(() => {
    frame.style.height = Math.ceil(document.body.getBoundingClientRect().height) + "px";
  }).observe(document.body);
}
