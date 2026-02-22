use std::f32::consts::PI;
use std::time::Instant;

use iced::gradient::Linear;
use iced::widget::canvas::{self, stroke, Cache, Canvas, Geometry, Path, Stroke};
use iced::window;
use iced::{mouse, Color};

use iced::{Element, Fill, Point, Rectangle, Renderer, Subscription, Theme};

pub fn main() -> iced::Result {
    iced::application(OSCPS::new, OSCPS::update, OSCPS::view)
        .subscription(OSCPS::subscription)
        .theme(Theme::Dark)
        .run()
}

struct OSCPS {
    start: Instant,
    cache: Cache,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    Tick,
}

impl OSCPS {
    fn new() -> Self {
        OSCPS {
            start: Instant::now(),
            cache: Cache::default(),
        }
    }

    fn update(&mut self, message: Message) {
        self.cache.clear();
    }

    fn view(&self) -> Element<'_, Message> {
        Canvas::new(self).width(Fill).height(Fill).into()
    }

    fn subscription(&self) -> Subscription<Message> {
        window::frames().map(|_| Message::Tick)
    }
}

impl Default for OSCPS {
    fn default() -> Self {
        OSCPS::new()
    }
}

impl<Message> canvas::Program<Message> for OSCPS {
    type State = ();

    fn draw(
        &self,
        _state: &Self::State,
        renderer: &Renderer,
        theme: &Theme,
        bounds: Rectangle,
        cursor: mouse::Cursor,
    ) -> Vec<Geometry> {
        let geometry = self.cache.draw(renderer, bounds.size(), |frame| {
            let palette = theme.palette();

            let center = frame.center();
            let radius = frame.width().min(frame.height()) / 5.0;

            let start = Point::new(center.x, center.y - radius); //  Top, stationary point

            // One full rotation every second
            let angle = (self.start.elapsed().as_millis() % 10_000) as f32 / 10_000.0 * 2.0 * PI;

            let end = match cursor {
                mouse::Cursor::Available(point) => point,
                mouse::Cursor::Levitating(point) => point,
                mouse::Cursor::Unavailable => Point::new(center.x, center.y),
            };

            // let end = Point::new(
            //     center.x + radius * angle.cos(),
            //     center.y + radius * angle.sin(),
            // );

            // Draw the end point circles
            let circles = Path::new(|b| {
                b.circle(start, 10.0);
                b.move_to(end);
                b.circle(end, 10.0);
            });

            let path = Path::new(|b| {
                b.move_to(start);
                b.arc_to(center, end, 50.0);
                b.line_to(end);
            });

            frame.stroke(
                &path,
                Stroke {
                    // style: stroke::Style::Solid(palette.text),
                    style: stroke::Style::Solid(Color::BLACK),
                    width: 10.0,
                    ..Stroke::default()
                },
            );

            frame.fill(&circles, palette.danger);
        });

        vec![geometry]
    }
}
