"use strict";

const fs = require("fs");
const vm = require("vm");

function buildHarness(options = {}) {
  const handlers = {};
  const shinyEvents = [];
  const copied = [];
  let fallbackCalls = 0;

  const document = {
    readyState: "loading",
    hidden: false,
    body: {
      appendChild: () => {},
      removeChild: () => {}
    },
    addEventListener: (name, callback) => { handlers[name] = callback; },
    getElementById: (id) => id === "chron-out_age_depth"
      ? { textContent: options.text === undefined ? "Plot()\n{\n};" : options.text }
      : null,
    createElement: () => ({
      value: "",
      style: {},
      setAttribute: () => {},
      select: () => {},
      setSelectionRange: () => {}
    }),
    execCommand: (command) => {
      if (command !== "copy") throw new Error("Unexpected command");
      fallbackCalls += 1;
      return options.fallbackSucceeds !== false;
    }
  };

  const context = {
    Promise,
    Date,
    window: {
      location: { pathname: "/" },
      top: {},
      setInterval: () => 1,
      isSecureContext: options.secure !== false,
      navigator: options.modern === false ? {} : {
        clipboard: {
          writeText: (text) => {
            copied.push(text);
            return options.modernRejects
              ? Promise.reject(new Error("denied"))
              : Promise.resolve();
          }
        }
      },
      Shiny: {
        setInputValue: (name, value) => shinyEvents.push({ name, value })
      }
    },
    document
  };

  vm.runInNewContext(
    fs.readFileSync("www/custom.js", "utf8"),
    context,
    { filename: "www/custom.js" }
  );

  return {
    click: async () => {
      handlers.click({
        target: {
          closest: (selector) => selector === ".varg-copy-to-clipboard"
            ? {
                getAttribute: (name) => name === "data-copy-target"
                  ? "chron-out_age_depth"
                  : "chron-copy_age_depth_status"
              }
            : null
        }
      });
      await new Promise((resolve) => setImmediate(resolve));
    },
    copied,
    shinyEvents,
    fallbackCalls: () => fallbackCalls
  };
}

(async function () {
  const modern = buildHarness();
  await modern.click();
  if (modern.copied[0] !== "Plot()\n{\n};") {
    throw new Error("Modern Clipboard API did not receive generated code");
  }
  if (!modern.shinyEvents.at(-1).value.ok ||
      modern.shinyEvents.at(-1).value.method !== "clipboard") {
    throw new Error("Modern copy success was not reported to Shiny");
  }

  const fallback = buildHarness({ modernRejects: true });
  await fallback.click();
  if (fallback.fallbackCalls() !== 1 ||
      fallback.shinyEvents.at(-1).value.method !== "fallback") {
    throw new Error("Rejected Clipboard API did not use the copy fallback");
  }

  const empty = buildHarness({ text: "  " });
  await empty.click();
  if (empty.copied.length !== 0 ||
      empty.shinyEvents.at(-1).value.reason !== "empty") {
    throw new Error("Empty generated code was not rejected before copying");
  }

  console.log("chronology_clipboard_ok");
})();
