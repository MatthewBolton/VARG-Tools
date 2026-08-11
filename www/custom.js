(function () {
  "use strict";

  var idleWarningMs = 50 * 60 * 1000;
  var idleTimeoutMs = 60 * 60 * 1000;
  var idleCheckMs = 15 * 1000;
  var lastActivityAt = Date.now();
  var warningSent = false;
  var timeoutSent = false;
  var heavyJobActive = false;
  var heavyJobHandlerRegistered = false;

  function isHostedSession() {
    var pathParts = window.location.pathname.split("/");
    return pathParts[1] === "app_proxy" && Boolean(pathParts[2]);
  }

  function sendShinyEvent(name) {
    if (!window.Shiny || typeof window.Shiny.setInputValue !== "function") {
      return false;
    }
    window.Shiny.setInputValue(name, Date.now(), { priority: "event" });
    return true;
  }

  function recordActivity() {
    if (timeoutSent) {
      return;
    }

    lastActivityAt = Date.now();
    if (warningSent) {
      warningSent = false;
      sendShinyEvent("session_idle_resumed");
    }
  }

  function setHeavyJobState(message) {
    var wasActive = heavyJobActive;
    heavyJobActive = Boolean(message && message.active);

    if (heavyJobActive || wasActive) {
      recordActivity();
    }
  }

  function registerHeavyJobHandler() {
    if (heavyJobHandlerRegistered || !window.Shiny ||
        typeof window.Shiny.addCustomMessageHandler !== "function") {
      return;
    }

    window.Shiny.addCustomMessageHandler("sessionHeavyJobState", setHeavyJobState);
    heavyJobHandlerRegistered = true;
  }

  function checkIdleTime() {
    if (timeoutSent) {
      return;
    }

    if (heavyJobActive) {
      lastActivityAt = Date.now();
      return;
    }

    var idleForMs = Date.now() - lastActivityAt;
    if (idleForMs >= idleTimeoutMs) {
      if (sendShinyEvent("session_idle_timeout")) {
        timeoutSent = true;
      }
      return;
    }

    if (idleForMs >= idleWarningMs && !warningSent) {
      if (sendShinyEvent("session_idle_warning")) {
        warningSent = true;
      }
    }
  }

  function startIdleTracking() {
    if (!isHostedSession()) {
      return;
    }

    registerHeavyJobHandler();
    document.addEventListener("shiny:connected", registerHeavyJobHandler);

    ["pointerdown", "pointermove", "keydown", "touchstart", "wheel", "scroll"].forEach(function (eventName) {
      document.addEventListener(eventName, recordActivity, { passive: true });
    });
    document.addEventListener("visibilitychange", function () {
      if (!document.hidden) {
        recordActivity();
      }
    });

    window.setInterval(checkIdleTime, idleCheckMs);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", startIdleTracking);
  } else {
    startIdleTracking();
  }
})();

(function () {
  "use strict";

  window.VARGTools = window.VARGTools || {};

  window.VARGTools.stopHostedSession = function (hostWindow) {
    var host = hostWindow || window.top;
    var stopHandler = host && host.Shiny && host.Shiny.instances &&
      host.Shiny.instances.eventHandlers &&
      host.Shiny.instances.eventHandlers.onStopApp;

    if (typeof stopHandler === "function") {
      var originalConfirm = host.confirm;
      host.confirm = function () { return true; };
      try {
        return Promise.resolve(stopHandler.call(host.Shiny.instances.eventHandlers));
      } catch (error) {
        return Promise.reject(error);
      } finally {
        host.confirm = originalConfirm;
      }
    }

    var stopControl = host && host.document &&
      host.document.querySelector(".app-stop-current-btn");
    if (stopControl) {
      var fallbackConfirm = host.confirm;
      host.confirm = function () { return true; };
      try {
        stopControl.click();
        return Promise.resolve();
      } catch (error) {
        return Promise.reject(error);
      } finally {
        host.confirm = fallbackConfirm;
      }
    }

    return Promise.reject(new Error("ShinyProxy stop handler was not found"));
  };
})();
